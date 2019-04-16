import lib.componentdetection.pcore as pc
from pyspark.sql import SparkSession
import sys
import shutil
import os
import time
from numpy import *

def create_spark_Session():
    # create spark session
    spark = SparkSession \
        .builder \
        .appName("ParsecComponentDetection") \
        .getOrCreate()
    return spark

def format_peak_output(p):
    # s = sorted(x, key=lambda p: p.Intensity, reverse=True)
    return "%f,%f,%f,%s,%i,%f" % \
           (p.RT, p.MZ, 0.3, "+MS", 1, p.Intensity)

def format_trace_output(t):
    # s = sorted(x, key=lambda p: p.Intensity, reverse=True)
    return "%s" % \
           (t.HighestIntensityMass)

def format_feature_output(f):
    # s = sorted(x, key=lambda p: p.Intensity, reverse=True)
    return "%f,%f,%f,%s,%i,%f,%f,%s" % \
           (f.RT, f.A0MZ, 0.3, "+MS", 1, f.Intensity, f.RtOffset,f.FileID)

def run(spark,i,o,sn,err,ref):

    #To connect to data in S3 with AWS creds
    #sc = spark.sparkContext
    #sc._jsc.hadoopConfiguration().set("fs.s3n.awsAccessKeyId", "ACCESSKEYID HERE")
    #sc._jsc.hadoopConfiguration().set("fs.s3n.awsSecretAccessKey", "ACCESSKEY HERE")

    output_dir_exists = os.path.exists(o)
    if output_dir_exists:
        print("removing the existing directory : ", o)
        shutil.rmtree(o)
    start = time.time()
    # Detect all features from all files
    features = spark.read.parquet(i) \
        .rdd.map(lambda scan: pc.kv_scan(scan)) \
        .filter(lambda scan: scan[1].MsOrder == 1).groupByKey() \
        .flatMap(lambda file: pc.bin_masses(file)) \
        .flatMap(lambda fileMassBin: pc.create_mass_traces(fileMassBin)) \
        .flatMap(lambda trace: pc.detect_peaks(trace[1], sn)) \
        .map(lambda p: pc.map_peak(p,0.1)) \
        .reduceByKey(lambda p1, p2: p1 + p2) \
        .map(lambda p: pc.detect_smallmolecule_features(p,err,False)) \
        .flatMap(lambda f: pc.flatten_features(f))

    features.map(lambda f: format_feature_output(f)) \
        .saveAsTextFile(o)

    end = round(time.time(),2)
    print("Parsec Workflow: ", end - start, "sec")
    count_results(spark,o)

def count_results(spark,input):
    print("# of Objects ",spark.read.text(input).count())

if __name__ == "__main__":

    s = create_spark_Session()
    run(s,sys.argv[1],sys.argv[2],1,50,'')

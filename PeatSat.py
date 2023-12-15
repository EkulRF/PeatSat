from PreProcessing.SAR_retrieval import *

bbox = [-4.4375, 56.0702778, -3.7752778, 58.3775]
time_range = ["2020-12-01","2022-12-31"]
sitename = 'FlowCountry'

S1_imagery = SAR_retrieval(bbox, time_range, sitename)


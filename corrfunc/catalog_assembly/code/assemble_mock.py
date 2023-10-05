from functions import assemble_data, assemble_rand
import sys

mock_id = sys.argv[1]

print("Starting Mock %s" % mock_id)
assemble_data("../data/EZmock_eBOSS_LRG_NGC_v7_%s_Nz-matched_mike.fits" % mock_id,
                         "../output/EZmock_eBOSS_LRG_NGC_v7_%s_Nz-matched_mike.dat" % mock_id, mock=True)
assemble_data("../data/EZmock_eBOSS_LRG_SGC_v7_%s_Nz-matched_mike.fits" % mock_id,
                         "../output/EZmock_eBOSS_LRG_SGC_v7_%s_Nz-matched_mike.dat" % mock_id, mock=True)
assemble_rand("../data/EZmock_eboss_LRG_NGC_mock-%s_rand.fits" % mock_id,
                         "../output/EZmock_eBOSS_LRG_NGC_v7_%s_Nz-matched_mike.rand" % mock_id, mock=True)
assemble_rand("../data/EZmock_eboss_LRG_SGC_mock-%s_rand.fits" % mock_id,
                         "../output/EZmock_eBOSS_LRG_SGC_v7_%s_Nz-matched_mike.rand" % mock_id, mock=True)


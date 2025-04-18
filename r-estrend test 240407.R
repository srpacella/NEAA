# Test using USGS R-Estrend package for seasonal kendall analyses
# This does not seem to work - Github notes it is currently an orphaned package
# https://code.usgs.gov/water/analysis-tools/restrend


remotes::install_gitlab("water/analysis-tools/smwrData",
                        host = "code.usgs.gov")
remotes::install_gitlab("water/analysis-tools/smwrBase",
                        host = "code.usgs.gov")
remotes::install_gitlab("water/analysis-tools/smwrGraphs",
                        host = "code.usgs.gov")
remotes::install_gitlab("water/analysis-tools/smwrStats",
                        host = "code.usgs.gov")
remotes::install_gitlab("water/analysis-tools/smwrQW",
                        host = "code.usgs.gov")
remotes::install_gitlab("water/analysis-tools/restrend",
                        host = "code.usgs.gov")

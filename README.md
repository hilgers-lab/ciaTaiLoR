# ciaTaiLoR
## R package for feature extraction and long-read processing
<!-- badges: start -->

  ![GitHub release (latest SemVer)](https://img.shields.io/github/v1.0.0/release/hilgers-lab/ciaTaiLoR)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-brightgreen)](https://github.com/hilgers-lab/ciaTaiLoR/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
  [![Downloads](https://img.shields.io/github/downloads/hilgers-lab/ciaTaiLoR/total)]()
  ![GitHub](https://img.shields.io/github/license/hilgers-lab/ciaTaiLoR)
  <!-- badges: end -->

    LASER determines the regulatory connections between exons, 5' ends, and 3' ends by analyzing every read as a complete transcript and using multinomial testing to evaluate the frequency of co-occurrence among these features.


  ### Installation

  ```
  install.packages("devtools")
  devtools::install_github("hilgers-lab/ciaTaiLoR")
  ```

  ```
  library(ciaTaiLoR)
  vignette("ciaTaiLoR")
  ```
  ### 3'end correction and long-read processing

  For 3'end correction and long-read processing and example vignette is available [here](https://hilgers-lab.github.io/polyADataBase/docs/polyADatabase.html)
# Release

Initial Release 0.1.0

Release date: 10th March 2023
This release corresponds to all associated processing functions and transcriptome used in by Alfonso-Gonzalez et al. manuscript

## Contact

Developer Carlos Alfonso-Gonzalez. For questions or feedback you can contact:

alfonso@ie-freiburg.mpg.de

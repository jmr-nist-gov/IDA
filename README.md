# Isotope Dilution Assistant

_current version 1.0; last updated 20220419_
_previous version 0.9; last updated 20180202_

Preprocessing of data in a consistent, curated manner is a key step toward any data product. Friction occurs often due to the need for balance between automation and expert assessment. One measurement technique that has suffered recently from such friction is time-resolved isotope dilution mass spectrometry (TR-IDMS). This is an analytical method in which the internal standard is an isotope of the same element. It has been classified as a primary higher-order method because of its inherent accuracy and is the method of choice for National Metrology Institutes (NMI) in value assigning Certified Reference Materials (CRM) for trace elements. Quantitation can be accomplished simply by measurement of the equilibrated isotope ratios using a mass spectrometer. The use of IDMS for value assignments at NIST has been declining recently because it requires specialist scientific expertise and involves extensive data processing and calculations, which can be quite labor intensive.

This data processing tool was constructed to directly address the time-consuming process of manual manipulation of text files into a calculation spreadsheet; "IDA" automates the process to rapidly calculate the mean, standard deviation, and relative standard deviation (RSD) of isotopic ratios.  It reads the basic comma separated value file generated by [SOFTWARE] for isotopic measurements of inorganic analytes.  IDA is a [Shiny](https://shiny.rstudio.com/) application, currently intended to be run from a local copy of RStudio.  The interface allows for changing of parameters to identify the signal stability region and manual override by brushing directly on the signal ratio plots.  Additionally, it flags potentially troublesome samples for the user's attention.  For best results and to use all features, [clone](https://gitlab.nist.gov/gitlab/jmr3/IDA.git) or [download](https://gitlab.nist.gov/gitlab/jmr3/IDA/-/archive/master/IDA-master.zip) from GitLab to a local instance and launch directly from RStudio.

`git clone https://gitlab.nist.gov/gitlab/jmr3/IDA.git`

Open the project "IDA" in RStudio. Open one of `global.R`, `ui.R`, or `server.R` within the `inst` directory. If your version of RStudio supports launching Shiny apps, there will be a "Run App" button toward the top right of the panel containing the opened file. IDA will verify you have the packages necessary to run. Full instructions are provided in the app, or see the [instructions page](instructions.html). An [example data file](inst/SRM 2778.csv) is provided that can also be loaded directly in the application.

Updates in Version 1.0 include bug fixes, better feedback, and improved handling of archival data sets.

IDA automates the entire process and provides downloads of result data including plots to support report writing. There is limited support for data sets containing unexpected data (e.g. probe failures resulting in data with unanticipated properties).



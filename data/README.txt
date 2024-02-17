In general, this folder is holding Julia serialization objects that are use for hardcoded calibrations of the detectors for APOGEE. 
Calibrations of this form are usually either done per telescope, chip, and/or fiber and so are stored as nested dictionaries with keys in that order,
down to the level of granualarity required. These files are small and git versioned, so we hope this exposes and makes tracking changes in hardcoded
calibrations in the pipeline easier. 
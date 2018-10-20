# Data Files

### Bethel Test Fishery CPUE Data

**File Name**: `BTF_data.txt`

**Source**: Data collected by the Alaska Department of Fish and Game, available [here](http://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareakuskokwim.btf).

**Variables**:

*  `year`: year of sampling
*  `day`: day of sampling (1 = 6/1 always, last is 85 = 8/24 always)
*  `doy`: day of year (varies depending on leap vs. non-leap year)
*  `date`: M/DD formated date
*  `cpue`: daily CPUE value
*  `ccpue`: cumulative CPUE value
*  `p.ccpue`: the fraction of end-of-season CCPUE that was caught by each date

---

### Total Estimated Abundance

**File Name**: `Total Run_Data.txt`

**Source**: Liller et al. [2018](http://www.adfg.alaska.gov/FedAidPDFs/RIR.3A.2018.04.pdf)

**Variables**:

*  `year`: the year of the estimate
*  `N`: the estimated drainage-wide Chinook salmon abundance in the Kuskokwim River.

---

### Observed and Forecast Median Run Date

**File Name**: `Total Run_Data.txt`

**Source**: Staton et al. [2017](https://www.sciencedirect.com/science/article/pii/S0165783617301248)

**Variables**:

*  `year`: the year of the timing estimates
*  `d50`: the doy of 50% of the run complete (observed)
*  `fcst_d50`: the mean run timing forecast that year
*  `fcst_se_d50`: the standard error of prediction for the run timing forecast that year

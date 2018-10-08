This is the online supplement for the article by B. Staton and M. Catalano published in CJFAS entitled:

"Bayesian information updating procedures for Pacific salmon run size indicators: Evaluation in the presence and absence of auxiliary migration timing information"

This supplement contains all of the code and data used in the analysis. 

Code Files:

NOTE: The online web portal to upload these files would not accept .R files so we have pasted the code into .txt files.
      To run the code, just copy and paste the code from the 3 *_Code.txt files into blank .R files and save them with the same name.

   Analysis.txt: contains the code that carries out the analysis

   Functions.txt: contains the functions we wrote to help carry out the analysis

   Plotting.txt: contains the code that makes all of the plots and tables shown in the article. 

Data Files (see article for full citations):

   2a_BTF_Data.txt: contains daily and cumulative CPUE observations from the Bethel Test Fishery (Bue and Lipka 2016). 
             The raw data may be accessed at: http://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareakuskokwim.btf

   2b_Cumulative Harvest_Data.txt: contains cumulative harvest by commercial and subsistence fisheries downstream of the Bethel Test Fishery.
                       Data were provided to B. Staton by N. Smith and Z. Liller (ADF&G) personally, 
		       but see Hamazaki 2011 for the approach for estimating subistence harvest

   2c_Regression_Data.csv: contains the vulnerable run size passing the Bethel test fishery and end-of-season CPUE from the Bethel test fishery in each year

   2d_Run Timing_Data.csv: contains the fitted and forecasted values of the run timing curve used by the in-season updating method
                          that included run timing forecast information. Methods for obtaining these quantities are described in Staton et al. (2017)

   2e_Total Run_Data.csv: contains the total estimated run size in each year (1976-2017) from the run reconstruction model, as presented in Smith and Liller (2018)

Running the Analysis:

   To run the analysis, open the Analysis.R file, and change line 24 to be the directory containing all of these files. You can then run the whole script.
   An output folder will be created, and the necessary output to make the plots in Plotting.R will be dumped there. The analysis should take ~40 minutes to complete.

   Investigation of how the approach works will be best conducted by diving in to the Functions.R file. A description of the meaning of each function argument is provided therein.

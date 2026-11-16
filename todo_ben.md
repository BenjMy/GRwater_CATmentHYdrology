- [ ] where is the ETp time serie input data located? and why its not defined in the config file?
- [ ] add a Ks with 2 layers scenarii for the real data (Agramon)  
- [ ] write a more generic config file where all the input data path are defined and import it from the pipeline... 
The idea is is to make it more generic to welcome new data that are not agramon catchements
- [ ] the scenario_1_plot3_outlet-closed_bottom-closed_spinup-on0c_pmin-m15_ic-wt_forcing-etpts-wet-fid2 shows wrong monthly LAI in veg_map_2016_monthly
- [ ] integrate plot_fig8_ETa_ETp_vs_psi into the pipeline or an existing plot file
- [ ] improve python run_paper_scenarios.py. I want it to test different hypotheses. Create a python file for each and then a run_all_hypotheses.py file that call each of them. 
  - hyp1_layers_Ks.py would test different layers configurations (1L, 2L, various Ks way handling etc ..)
  - hyp2_LAI_SCF.py would test different CATHY version (with/without LAI SCF and PMIN values -15/-150/-1e35) 
  - hyp3_ic.py would test different configureation of ic (depth of water table)
  - hyp4_SW.py would test different ETp inputs (from wet ShuterWallace to dry SW)
  - with/without spinup... etc 
  Each must create a single folder with all the run inside 
  refactors plotters accordindly so you plot interpret for each hyp to test. 
  ideally run_all_hypotheses.py must contains a dry-run call so I can visualise quickly the runs 
  ideally run_all_hypotheses.py must also call the plotters too 
  
  each files hyp* must contains the validation steps versus EB data (.netcdf). See for instance analysis_ET_compare.py
  
  update the readme accordindly 
  
  you can still based all the studies/scenarii on the .toml (thus improve/complete it)
  
- [ ] remove batch_run_20261116_144740_plot3_pmin-50--1e35_outlet-both_bottom-both_spinup-off_ic-both.log writting
- [ ] remove to_run_paper.md (deprecated)

- [ ] improve lib filenames nomenclature for readability.

#!/usr/bin/python
#
# Prepare templates for output files from offline radiative transfer calculations 
#   suitable for publishing on the Earth System Grid 
#
# Robert Pincus, Robert.Pincus@colorado.edu, 2016-2017
#
# ---------------------------------------------------------------------------------
from netCDF4 import Dataset  
# ---------------------------------------------------------------------------------
# Copy a variable and all its attributes from one netCDF file to another 
#
def copyVar(nc_in, nc_out, name, newname=None) :   
	if newname is None : 
		newname = name
	nc_out.createVariable(newname, nc_in.variables[name].dtype, nc_in.variables[name].dimensions) 
	nc_out.variables[newname].setncatts(nc_in.variables[name].__dict__)
	nc_out.variables[newname][:] = nc_in.variables[name][:] 
# ---------------------------------------------------------------------------------

atmos_file = Dataset('multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-0-4_none.nc', mode='r') 
# Available from https://www.earthsystemcog.org/projects/rfmip/resources/
# or from https://esgf-node.llnl.gov/search/input4mips/ ; search for "RFMIP" 

short_names = ['rlu','rsu', 'rld', 'rsd'] 
stand_names = ['upwelling_longwave_flux_in_air','upwelling_shortwave_flux_in_air', 
               'downwelling_longwave_flux_in_air','downwelling_shortwave_flux_in_air'] 

# Attributes are take from https://docs.google.com/document/d/1h0r8RZr_f3-8egBMMh7aqLwy3snpD6_MrDz1q8n5XUk/edit
# Data reference syntax attributes 
drs_attrs = { 
  "activity_id"  :"RFMIP",   # (from CMIP6_activity_id.json)
  "product"      :"model_output",
  "experiment_id":"rad-irf", # (from CMIP6_experiment_id.json)
  "table_id"     :"Efx",     # (per http://clipc-services.ceda.ac.uk/dreq/u/efc0de22-5629-11e6-9079-ac72891c3257.html) 
  "frequency"    :"fx", 
  "sub_experiment_id":"none"} 

expt_attrs = { 
  "Conventions"         :"CF-1.7 CMIP-6.0", 
  "experiment"          :"rad_irf", 
  "sub_experiment"      :"none", 
  "realization_index"   :1, 
  "initialization_index":1,
  "nominal_resolution"  :"50 km", 
  "grid"                :"columns sampled from ERA-Interim, radiative fluxes computed indepednently"} 
  
# Further required attributes, uniform across submissions 
std_attrs = { 
  "data_specs_version":"1.00.12", 
  "forcing_index":1, # This values follows page 2074 in https://dx.doi.org/10.5194/gmd-10-2057-2017
                     # 1 = calculations uses all available greenhouse gases
                     # 2 = calculation uses CO2, CH4, N2O, CFC11eq
                     # 3 = calculation uses CO2, CH4, N2O, CFC12eq, HFC-134eq
  "physics_index":1} 

# Model/institution specific attributes 
model_attrs = {
  "institution_id"  :"",
  "source_id"       :"",
  "further_info_url":"",
  "license"         :""} 

# Submission attrs 
sub_attrs = {
  "creation_date":"",
  "variant_label":"", 
  "version"      :""} 

for short, std in zip(short_names, stand_names) : 
    out_file_name = short + "_template.nc"
    print('Creating ' + out_file_name)
    out_file = Dataset(out_file_name, mode='w', FORMAT='NETCDF4_CLASSIC')
    out_file.setncatts(drs_attrs) 
    out_file.setncatts(std_attrs) 
    out_file.setncatts(expt_attrs) 
    out_file.setncatts(model_attrs) 
    out_file.setncatts(sub_attrs) 
    d = out_file.createDimension('expt',  atmos_file.dimensions['expt'].size) 
    d = out_file.createDimension('site',  atmos_file.dimensions['site'].size) 
    d = out_file.createDimension('level', atmos_file.dimensions['level'].size) 
    copyVar(atmos_file, out_file, 'lat') 
    copyVar(atmos_file, out_file, 'lon') 
    copyVar(atmos_file, out_file, 'time') 
    copyVar(atmos_file, out_file, 'pres_level', 'plev') 
    v = out_file.createVariable(short,   'f4', ('expt', 'site', 'level')) 
    v.setncatts({'variable_id'  :short,
                 'standard_name':std, 
                 'units'        :'W m-2', 
                 'coordinates'  : 'lon lat time'}) 
    copyVar(atmos_file, out_file, 'profile_weight') 
    out_file.close()



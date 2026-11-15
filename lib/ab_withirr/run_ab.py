"""A/B: rain+ETp combined in atmbc (A) vs rain via irrigation_map with atmbc=-ETp (B). Flat 16x16 case."""
import os, sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
mode = sys.argv[1]; out = sys.argv[2]
from pyCATHY import cathy_tools
from pyCATHY.irrigation import write_irrigation_map
NROW=NCOL=16; DX=DY=25.0; NSTR=12; BASE=6.0; DAY=86400.0; DAYS=int(os.environ.get("DAYS","60"))
ZR=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.10,0.12,0.18,0.24]
rng=np.random.default_rng(1)
rain_mm=np.zeros(DAYS); rain_mm[::10]=rng.uniform(15,30,len(rain_mm[::10]))   # storm every 10 d
et_mm=np.full(DAYS,5.0)*np.where(rain_mm>0,0.5,1.0)
t_end=DAYS*DAY; times=np.append(np.arange(DAYS)*DAY,t_end)
rain_ms=np.append(rain_mm,rain_mm[-1])/1000/DAY; et_ms=np.append(et_mm,et_mm[-1])/1000/DAY
simu=cathy_tools.CATHY(dirName=out,prj_name="ab",version="withIrr")
ii,jj=np.meshgrid(np.arange(NROW),np.arange(NCOL),indexing="ij")
DEM=100.0-1e-4*DX*(ii+jj)/2.0
simu.update_dem_parameters(delta_x=DX,delta_y=DY,nstr=NSTR,base=BASE,zratio=ZR)
simu.update_prepo_inputs(DEM,delta_x=DX,delta_y=DY)
simu.update_zone(np.ones((NROW,NCOL),dtype=int)); simu.update_veg_map(np.ones((NROW,NCOL))); simu.update_soil()
if mode=="A":
    simu.update_atmbc(HSPATM=1,IETO=1,time=times,VALUE=[rain_ms,et_ms])
    write_irrigation_map(simu._input_path("irrigation_map"),times[:1],np.zeros((1,NROW,NCOL)),irr_sw=1.0,irr_flag=0)
else:
    simu.update_atmbc(HSPATM=1,IETO=1,time=times,VALUE=[np.zeros_like(rain_ms),et_ms])
    rates=np.broadcast_to(rain_ms[:-1,None,None],(DAYS,NROW,NCOL)).copy()
    write_irrigation_map(simu._input_path("irrigation_map"),times[:-1],rates,irr_sw=1.01,irr_flag=1)
simu.update_ic(INDP=0,IPOND=0,pressure_head_ini=-5.0)
simu.update_parm(TIMPRTi=list(np.arange(1,DAYS+1)*DAY),VTKF=0,TMAX=t_end,IPRT=2)
simu.run_preprocessor(verbose=False)
simu.run_processor(IPRT1=2,TRAFLAG=0,DELTAT=60.0,DTMIN=1e-4,DTMAX=21600.0,TMAX=t_end,verbose=False)
np.savez(os.path.join(out,"forcing.npz"),rain_mm=rain_mm,et_mm=et_mm)
print("RUN DONE",mode)

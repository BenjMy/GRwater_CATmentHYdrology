# RUN
# -------------------------------------------------

KS_VALUES=(1e-4)

ZROOT_SETS=(
    "0.1 0.3 1.0 1.5 2.5"   # baseline
    "2.5 2.5 2.5 2.5 2.5"   # deeper rooting
)

for ks in "${KS_VALUES[@]}"; do
    for zf in "${ZROOT_SETS[@]}"; do
        ./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -1e34 -15\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${ks}" \
            --zroot-factors ${zf}
    done
done



KS_VALUES=(1e-4)

# Each row = bare, sparse, moderate, dense, very-dense multiplier on --zroot
ZROOT_SETS=(
    "0.1 0.3 1.0 1.5 2.5"   # baseline
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)

for ks in "${KS_VALUES[@]}"; do
    for zf in "${ZROOT_SETS[@]}"; do
        ./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 2 \
            --pmin -1e34 -15\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${ks}" \
            --zroot-factors ${zf}
    done
done


KS_VALUES=(1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)

./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -15 -1e34\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS}



KS_VALUES=(1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)

./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -150\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS}




KS_VALUES=(1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)
./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -1e34\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS} \
            --pcana 5


KS_VALUES=(1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)
./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -50\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS} \
            --pcana 5



KS_VALUES=(6e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)
./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -50\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS} \

Good one
------
KS_VALUES=(1e-4 6e-4 1e-3)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)
./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -50 -1e35\
            --outlet both --bottom both \
            --ic both \
            --ic-pressure-head -100 \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS} \




KS_VALUES=(6e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)
./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -50 -1e35\
            --outlet closed --bottom closed \
            --ic both \
            --ic-wtposition -1 \
            --ic-pressure-head -100 \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS} \
            --scf 1



KS_VALUES=(6e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)
./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -50 -1e35\
            --outlet closed --bottom closed \
            --ic pressure \
            --ic-pressure-head -100 \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS} \
            --scf 1






KS_VALUES=(1e-6 1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)

./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 1\
            --pmin -1e34 -15\
            --outlet closed --bottom closed \
            --ic both \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS}



KS_VALUES=(1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)

./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -1e34 \
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES[0]}" \
            --zroot-factors ${ZROOT_SETS[0]} \
            --scf 0.5



KS_VALUES=(1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)

./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -50 \
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES[0]}" \
            --zroot-factors ${ZROOT_SETS[0]} \
            --scf 0






KS_VALUES=(1e-4)
ZROOT_SETS=(
    "0.3 0.6 2.0 3.0 5.0"   # deeper rooting
)

./run_scenarii_withETp_withLAI_ALL.sh \
            --dem-plot 3 \
            --pmin -1e34\
            --outlet closed --bottom closed \
            --ic pressure \
            --permx "${KS_VALUES}" \
            --zroot-factors ${ZROOT_SETS} \
            --spinup on \
            --spinup-cycles 3 \







# Improvements 
# -------------------------------------------------

## Ponded and saturated states 
are distinguished on the basis of a threshold parameter, hmin [m], that represents the minimum water
depth before surface routing can occur. In other words, the model allows surface water to flow only when the ponding head is equal
or larger than the threshold; hmin takes into account surface roughness and microtopography, starting from zero for perfectly smooth
surfaces and becoming larger for increasingly rough surfaces.

## Ks and LAI Thomas et al. 2021




# Bugs 
# -------------------------------------------------


# In the paper you read, I think (ten years have passed already!) I interpolated simulated ETa over time to better compare with the measured ETa, that's probably why you don't see ETa going to 0 in the plots.

As for SCF, it was intended to be used either at 0 (PMIN threshold approach for calculating ETa) or 1 (Feddes approach). I never  really tested it with intermediate values, as I suspect it would lead to numerical errors. We typically assume that the major contribution to ET is transpiration (SCF=1). This is a reasonable assumption except for arid ecosystems, where the vegetation cover is typically sparse and there is a combination of bare soil and plants. This also depends on the spatial resolution of course.


# bugs for dem 2


  File "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision/Supervision/Xela_Carracedo_Practicas_2026_data/Codigos_HydroModel/Agramon_withLAI_withETp.py", line 2971, in <module>
    run_pipeline(args)
    ~~~~~~~~~~~~^^^^^^
  File "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision/Supervision/Xela_Carracedo_Practicas_2026_data/Codigos_HydroModel/Agramon_withLAI_withETp.py", line 2207, in run_pipeline
    print(f"  Mesh x : {float(ds_mesh.x.min()):.1f} – {float(ds_mesh.x.max()):.1f}")
                              ~~~~~~~~~~~~~^^
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/xarray/core/_aggregations.py", line 2899, in min
    return self.reduce(
           ~~~~~~~~~~~^
        duck_array_ops.min,
        ^^^^^^^^^^^^^^^^^^^
    ...<3 lines>...
        **kwargs,
        ^^^^^^^^^
    )
    ^
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/xarray/core/dataarray.py", line 3875, in reduce
    var = self.variable.reduce(func, dim, axis, keep_attrs, keepdims, **kwargs)
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/xarray/core/variable.py", line 1712, in reduce
    result = super().reduce(
        func=func, dim=dim, axis=axis, keepdims=keepdims, **kwargs
    )
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/xarray/namedarray/core.py", line 921, in reduce
    data = func(self.data, **kwargs)
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/xarray/core/duck_array_ops.py", line 499, in f
    return func(values, axis=axis, **kwargs)
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/xarray/core/nanops.py", line 71, in nanmin
    return nputils.nanmin(a, axis=axis)
           ~~~~~~~~~~~~~~^^^^^^^^^^^^^^
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/xarray/core/nputils.py", line 242, in f
    result = getattr(npmodule, name)(values, axis=axis, **kwargs)
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/numpy/lib/_nanfunctions_impl.py", line 363, in nanmin
    res = np.amin(a, axis=axis, out=out, **kwargs)
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/numpy/_core/fromnumeric.py", line 3316, in amin
    return _wrapreduction(a, np.minimum, 'min', axis, None, out,
                          keepdims=keepdims, initial=initial, where=where)
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/numpy/_core/fromnumeric.py", line 81, in _wrapreduction
    return reduction(axis=axis, out=out, **passkwargs)
  File "/home/z0272571a@campus.csic.es/micromamba/envs/centum/lib/python3.14/site-packages/numpy/_core/_methods.py", line 45, in _amin
    return umr_minimum(a, axis, None, out, keepdims, initial, where)
ValueError: zero-size array to reduction operation minimum which has no identity


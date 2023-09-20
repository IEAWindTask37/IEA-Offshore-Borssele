With PyWake, run this example as:

```
python pywake_ex.py ../IEA37_Borssele_Regular_System.yaml 
python pywake_ex.py ../IEA37_Borssele_Irregular_System.yaml 
python opt.py ../IEA37_Borssele_Regular_System.yaml
```

With Floris, run this example as:

```
python floris_ex.py ../IEA37_Borssele_Regular_System.yaml 
python floris_ex.py ../IEA37_Borssele_Irregular_System.yaml 
```

Differences in computed AEP between PyWake and Floris result from different
approaches for linear interpolation of the wind rose when calculating a
finer resolution for the frequency matrix.

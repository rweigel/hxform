# hxtransform
Heliophysics coordinate transform 

To compile hxtransform python package, execute
```bash  
f2py -c Geopack_08_dp.for T96_01.for -m hxtransform
```

To test if compilation was successful, execute
``` python -c "import hxtransform as hx; help(hx) ```

Which should produce something similar to 
```
<module 'hxtransform' from '/Users/user/hxtransform/hxtransform.cpython-38-darwin.so'> 
``` 

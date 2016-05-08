pro pyspec1directory, dir

  print, ''
  print, 'Running PYSPECFIT for all objects in : '+dir
  print, ''

  spawn, 'ls *.pyspec > tmp.list'
  readcol, 'tmp.list', files, f = 'A'
  nfiles = n_elements(files)
  ids = strarr(nfiles)
  for ii = 0, nfiles - 1 do $
     ids[ii] = strmid(files[ii], 0, )

  file_mkdir, 
  mpirun -np 12 python fitspec.py 1090_1 INNER 1.50600 0.05 1090_results/inner1


end

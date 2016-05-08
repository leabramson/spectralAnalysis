pro pyspec1ID, ID, PA, z

  ID = strcompress(string(ID), /rem)
  PA = strcompress(string(PA), /rem)
  z = strcompress(string(z), /rem)
  
  print, ''
  print, 'Running PYSPECFIT for ID '+ID+' PA '+Pa
  print, ''

  file_mkdir, ID+'_results'
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' INNER '+z+' 0.05 '+ID+'_results/inner'+PA
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' INTER '+z+' 0.05 '+ID+'_results/inter'+PA
  spawn, 'mpirun -np 12 python fitspec.py '+ID+'_'+PA+' OUTER '+z+' 0.05 '+ID+'_results/outer'+PA

end

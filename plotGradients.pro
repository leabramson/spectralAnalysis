function getResults, dir

  spawn, 'ls '+dir+'/inner.analyze > tmpInner.list'
  spawn, 'ls '+dir+'/inter.analyze > tmpInter.list'
  spawn, 'ls '+dir+'/outer.analyze > tmpOuter.list'

  
  
  RETURN, results
end

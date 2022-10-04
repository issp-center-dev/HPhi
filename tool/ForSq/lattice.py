def func_site(all_i,list_org):
  Lx      = list_org[0]
  Ly      = list_org[1]
  Lz      = list_org[2]
  orb_num = list_org[3]

  orb  = all_i%orb_num
  site = (all_i-orb)/orb_num
  x    = site%Lx
  tmp  = (site-x)/Lx
  y    = tmp%Ly
  z    = (tmp-y)/Ly 
  list_site = [int(x),int(y),int(z),int(orb)]
  return list_site


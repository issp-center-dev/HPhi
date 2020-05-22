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

def func_strans(list_trans,list_site,list_org):
  x_j   = int((list_site[0]+list_trans[0])%list_org[0])
  y_j   = int((list_site[1]+list_trans[1])%list_org[1])
  z_j   = int((list_site[2]+list_trans[2])%list_org[2])
  all_j = list_trans[3]+(x_j+y_j*list_org[0]+z_j*list_org[0]*list_org[1])*list_org[3]
  return all_j

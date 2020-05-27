def func_input(list_input):
  with open('input.txt') as f:
    data = f.read()
  data = data.split("\n")
  dict_lat = {}
  for i in data:
    tmp = i.split()
    if len(tmp)>0:
      for i in list_input:
        if i == tmp[0]:
          dict_lat[i] = tmp[1]
  return dict_lat



def func_readdef(name):
  with open(name) as f:
    data = f.read()
  data = data.split("\n")
  tmp = data[1].split()
  print(tmp[0],tmp[1])
  return int(tmp[1])


def func_count(file_name):
    #[s] file name
    #print("  name of input file = ",file_name)
    #[e] file name
    with open(file_name) as f:
        data      = f.read()
        data      = data.split("\n")
        #print(len(data))
        #[s] count not empty elements
    cnt = 0
    for i in range(0,len(data)):
        if data[i]: # if data[i] is not empty
           cnt += 1
        #print(cnt)
    return cnt

def func_readpair(file_name,siteI,siteJ,intT1,intT2,para):
    #[s] file name
    print("  name of input file = ",file_name)
    #[e] file name
    with open(file_name) as f:
        data      = f.read()
        data      = data.split("\n")
        #print(len(data))
        #[s] count not empty elements
    cnt   = 0
    for i in range(0,len(data)):
        if data[i]: # if data[i] is not empty
           tmp        = data[i].split()
           siteI[cnt] = int(tmp[0])
           siteJ[cnt] = int(tmp[1])
           intT1[cnt]  =     tmp[2]
           intT2[cnt]  =     tmp[3]
           para[cnt]  = float(tmp[4])
           cnt       += 1
        #print(cnt)

def func_input(list_lat):
    with open('input.txt') as f:
        data = f.read()
    data = data.split("\n")
    #print('check input',data)
    #print(len(data))
    dict_lat = {}
  
    for i in data:
        tmp = i.split()
        if len(tmp)>0:
            for i in list_lat:
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

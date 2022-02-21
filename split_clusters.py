!#usr/bin/env/python3
# this script splits the zscore general table in different files that contain all the genes corresponding to each cluster #

# zscore table #
file1 = open("/home/arnau/tpm_zscore.csv","r")
# metadata with identifier+nºcluster #
file2 = open("/home/arnau/final_clusters_signed.txt","r")

# empty lists to load the input and output data #
lista_r=[]
lista_c=[]
lista_completo=[]
   
# splitter function #
def create_files():
    i=0
    for line in file1:
        # if i==0:
        #     break
    
        data=line.split('\t')
        iden=data.pop(0)
        registro=Registro(iden,data)
        lista_r.append(registro)
        i=i+1
    
    i=0    
    for line in file2:
        data=line.split('\t')
        cluster=Cluster(data[0],data[1])
        lista_c.append(cluster)
        i=i+1
    
    for registro in lista_r:
        cluster_text=''
        for cluster in lista_c:
            if cluster.id==registro.id:
                cluster_text=cluster.cluster
        completo=Completo(registro.id,registro.datos,cluster_text)        
        lista_completo.append(completo)
    
    for completo in lista_completo:
        
        datos_formated = '\t'.join([str(elem) for elem in completo.datos])
        
        with open('./{0}.csv'.format(completo.cluster),'a')  as f:
             f.write('{0}\t{1}'.format(completo.id,datos_formated))
       
        
# input and output objects #        
        
class Completo:
    def __init__(self, id, datos, cluster):
        self.id = id
        self.datos = datos
        self.cluster = cluster

class Cluster:
    def __init__(self, id, cluster):
        self.id = id
        self.cluster = cluster


class Registro:
    def __init__(self, id, datos):
        self.id = id
        self.datos = datos

if __name__ == "__main__":
    create_files()


  
    


# coding: utf-8

#    # **_2D Truss System Analysis_**

# In[128]:


#Modüllerin Sisteme Yüklenmesi.
import math
import numpy as np
import matplotlib.pyplot as plt


# In[129]:


#Veri Girişi
L1=4
A=737
F_kn=2
E_gpa=20


# In[130]:


#Kabuller
DOF=2


# In[131]:


#Dönüşümler
E=E_gpa*10**3   #N/mm^2
l=L1*1000     #mm
F=F_kn*1000    #N


# In[132]:


#Çerçeve elemanların YEREL rijitlik matrisini oluşturan fonksiyon
def eleman_yerel_rijitlik_matrisi(x1,y1,x2,y2,E,A):
        L=math.sqrt(((x2-x1)**2)+((y2-y1)**2))
        c=float((x2-x1)/L)
        s=float((y2-y1)/L)
        K_eleman_yerel_rijitlik_matrisi=(E*A/L)*(np.array([[c**2,c*s,-c**2,-c*s],
                                                           [c*s,s**2,-c*s,-s**2],
                                                           [-c**2,-c*s,c**2,c*s],
                                                           [-c*s,-s**2,c*s,s**2]]))
        return(K_eleman_yerel_rijitlik_matrisi)


# In[133]:


#Eleman İsimleri, Noktalar ve Koordinatların oluşması
frames={1:[1,2],2:[1,3],3:[2,3],4:[2,4],5:[3,4],6:[3,5],7:[4,5],8:[4,6],9:[5,6],10:[5,7],11:[6,7]}
nodes={1:[0,0],2:[(l/2),l],3:[l,0],4:[((3*l)/2),l],5:[(2*l),0],6:[((5*l)/2),l],7:[(3*l),0]}


# In[134]:


#Sistem Geometrisinin Notasyona Göre Çizilmesi
for item in frames.keys():
    plt.plot([(nodes[frames[item][0]][0]),(nodes[frames[item][1]][0])],[(nodes[frames[item][0]][1]),(nodes[frames[item][1]][1])],"d-")
    for nokta in nodes.keys():
        plt.annotate(nokta, xy=(nodes[nokta][0], nodes[nokta][1]), xytext=(nodes[nokta][0], nodes[nokta][1]),
            bbox=dict(boxstyle="circle",alpha=0.5,fc="yellow"),
            )
for el_adı in frames.keys():
    x=((nodes[frames[el_adı][1]][0])+(nodes[frames[el_adı][0]][0]))/2
    y=((nodes[frames[el_adı][1]][1])+(nodes[frames[el_adı][0]][1]))/2
    plt.annotate(el_adı, xy=(x,y), xytext=(x,y),
        bbox=dict(boxstyle="round",alpha=1,fc="red"),
        )
plt.title('Sistem Geometrisi')
plt.show()


# In[154]:


#Elemanların Yerel RijitLik matrislerinin Hesaplanması   
K_eleman_genel_rijitlik_matrisi={}
for item in frames.keys():
    K_eleman_genel_rijitlik_matrisi[item]=eleman_yerel_rijitlik_matrisi((nodes[frames[item][0]][0]),(nodes[frames[item][0]][1]),(nodes[frames[item][1]][0]),(nodes[frames[item][1]][1]),E,A)
    #print("\n",20*"=","{}. Elemanın yerel rijitlik matrisi:".format(item),20*"=","\n")
    #print("{}. elemanın yerel rijitlik matrisi: \n".format(item))
    #print(np.around(K_eleman_genel_rijitlik_matrisi[item],3))


# In[136]:


# GENEL RİJİTLİK MATRİSİNİN BOYUTLARININ HESABI   
K_genel_rijitlik_matrisi_boyut = DOF * len(nodes)
K_genel_rijitlik_matrisi = np.zeros((K_genel_rijitlik_matrisi_boyut,K_genel_rijitlik_matrisi_boyut))
#print(K_genel_rijitlik_matrisi)


# In[137]:


#GENEL RİJİTLİK MATRİSİNİ HARİTALANDIRIP OLUŞTURAN FONKSİYON    
def sistem_genel_rijitlik_matrisinin_oluşturulması (K_eleman_genel_rijitlik_matrisi,nodal_names,DOF,K_genel_rijitlik_matrisi_boyut):
    K_sanal = np.zeros((K_genel_rijitlik_matrisi_boyut,K_genel_rijitlik_matrisi_boyut))
    # i,i 
    index = nodal_names[0] * DOF - DOF 
    K_sanal[ index : index + DOF , index : index + DOF  ] =  K_eleman_genel_rijitlik_matrisi[ 0 : DOF  , 0 : DOF ]
    # i,j
    index_2 = nodal_names[1] * DOF - DOF 
    K_sanal[ index : index + DOF , index_2 : index_2 + DOF ] =  K_eleman_genel_rijitlik_matrisi[ 0 : DOF  , DOF : DOF + DOF ]
    # j,i
    K_sanal[ index_2 : index_2 + DOF , index : index + DOF ] =  K_eleman_genel_rijitlik_matrisi[ DOF : DOF + DOF  , 0 : DOF ]
    # j,j
    index = nodal_names[1] * DOF - DOF 
    K_sanal[ index : index + DOF , index : index + DOF  ] =  K_eleman_genel_rijitlik_matrisi[ DOF : DOF + DOF , DOF : DOF + DOF ]
    return( K_sanal )


# In[138]:


# GENEL RİJİTLİK MATRİSİNİN HESAPLANMASI
for item in frames.keys():
    K_genel_rijitlik_matrisi += sistem_genel_rijitlik_matrisinin_oluşturulması( K_eleman_genel_rijitlik_matrisi[item] , frames[item], DOF , K_genel_rijitlik_matrisi_boyut )
    #print( K_genel_rijitlik_matrisi )
#print("\n",20*"=","SİSTEMİN GENEL RİJİTLİK MATRİSİ",20*"=","\n")
#print( np.round(K_genel_rijitlik_matrisi,3))
#print("\n"+"="*50)


# In[139]:


#SINIR KOŞULLAR
sınır_koşullar=np.array([0,0,1,1,1,1,1,1,1,1,1,1,0,0])
K_genel_rijitlik_matrisi_sınır_koşul=K_genel_rijitlik_matrisi.copy()
for item,value in enumerate(sınır_koşullar):
    if value==0:
        K_genel_rijitlik_matrisi_sınır_koşul[item,:]=0
        K_genel_rijitlik_matrisi_sınır_koşul[:,item]=0
        K_genel_rijitlik_matrisi_sınır_koşul[item,item]=1
#print("\n",20*"=","SİSTEMİN GENEL İNDİRGENMİŞ RİJİTLİK MATRİSİ",20*"=","\n")
#print(np.around(K_genel_rijitlik_matrisi_sınır_koşul))
#print("\n"+"="*50)
    


# In[140]:


#DEPLASMAN HESABI İÇİN Po IN OLUŞTURULMASI
Po=np.array([[0],[0],[0],[-F],[0],[0],[0],[-2*F],[0],[0],[0],[-F],[0],[0]])


# In[141]:


#ANALİZ İÇİN SINIR KOŞULLARA GÖRE İNDİRGENMİŞ Po VEKTÖRÜ
Po_sınır_koşul=Po.copy()
for item,value in enumerate(sınır_koşullar):
    if value==0:
        Po_sınır_koşul[item]=0


# In[142]:


#SİSTEM DEPLASMAN VEKTÖRÜNÜN HESAPLANMASI
U_deplasman_vektörü=np.linalg.solve(K_genel_rijitlik_matrisi_sınır_koşul,Po_sınır_koşul)
#print("SİSTEMİN FİKTİF YATAY YÜKLER ALTINDA DEPLASMANI")
print(np.round(U_deplasman_vektörü,4))
#print("\n"+"="*50)


# In[143]:


#Sistemin Mesnet Reaksiyonları ve Üzerinde Bulunan Yük Vektörleri
P_sistem=K_genel_rijitlik_matrisi.dot(U_deplasman_vektörü)
print(np.around(P_sistem,3))


# In[144]:


#Sistemin Genel Deplasman Vektörünün eleman genel vektörüne indirgenmesi
deplasman_eleman_vektör_sözlüğü_genel={}
for item in frames.keys():
    vektör=np.zeros([4,1])
    i=frames[item][0]
    j=frames[item][1]
    vektör[0]=U_deplasman_vektörü[(i*DOF)-2]
    vektör[1]=U_deplasman_vektörü[(i*DOF)-1]
    vektör[2]=U_deplasman_vektörü[(j*DOF)-2]
    vektör[3]=U_deplasman_vektörü[(j*DOF)-1]
    deplasman_eleman_vektör_sözlüğü_genel[item]=vektör


# In[145]:


#eleman genel deplasman vektörünün, eleman yerel deplasman vektörüne indirgenmesi
deplasman_vektörü_eleman_yerel={}
for item in frames.keys():
    deplasman_vektörü_eleman_yerel[item]=deplasman_eleman_vektör_sözlüğü_genel[item]


# In[146]:


#ELeman yerel yük vektörünün oluşturulması
eleman_yerel_yük_vektörü={}
for item in frames.keys():
    eleman_yerel_yük_vektörü[item]=(K_eleman_genel_rijitlik_matrisi[item].dot(deplasman_vektörü_eleman_yerel[item]))


# In[147]:


#ELeman yerel yük vektörünün oluşturulması
eleman_yerel_yük_vektörü={}
for item in frames.keys():
    aşama=(K_eleman_genel_rijitlik_matrisi[item].dot(deplasman_vektörü_eleman_yerel[item]))   
    eleman_yerel_yük_vektörü[item]=(aşama)
    #print("\n",40*"-")
    #print("\n{}. Elemanın Uç Kuvvetlerinin Vektörü (birim: Newton , N.mm):".format(item),"\n")
    #print(np.around(eleman_yerel_yük_vektörü[item],2))
    

    


# In[148]:


for item in frames.keys():
    x=eleman_yerel_yük_vektörü[item][0]
    if x>0:
        işaret=("-")
    else:
        işaret=("+")
    P_iç=((eleman_yerel_yük_vektörü[item][0])**2+(eleman_yerel_yük_vektörü[item][1])**2)**0.5
    print("\n",40*"-")
    print("{}.elemanın eksenel kuvveti:".format(item))
    print(str(işaret),float(np.around(P_iç,2)))


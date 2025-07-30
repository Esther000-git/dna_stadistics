# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 18:58:39 2022

@author: esther bermejo martinez
"""
from tkinter import Frame,Label,Tk,Entry,StringVar,ttk,DoubleVar,Scrollbar,Text,BOTTOM,E,W,S,N,END

from tkinter import Button as B


raiz=Tk()
raiz.title("DNA Stadistics")
raiz.iconbitmap("C:\\Users\\esthe\\Documentos_Proyecto ADN\\DNA_logo.ico")
raiz.config(bg="white")
raiz.geometry('1200x600')

frame=Frame()


p_GC=DoubleVar()    # Variables porcentajes. 
p_AT=DoubleVar()
p_AU=DoubleVar()
p_U=DoubleVar()
p_C=DoubleVar()
p_G=DoubleVar()
p_A=DoubleVar()
p_T=DoubleVar()

entry_PROTEINA=StringVar()   # Variables texto
entry_cadSilvestre=StringVar()
resultados_salida=StringVar()


raiz.columnconfigure(0, weight=0)
raiz.columnconfigure(1, weight=1)
raiz.columnconfigure(2,weight=1)
raiz.columnconfigure(3,weight=1)
raiz.columnconfigure(4,weight=1)
raiz.columnconfigure(5,weight=1)
raiz.columnconfigure(6,weight=1)
raiz.columnconfigure(7,weight=1)
raiz.columnconfigure(8,weight=1)
raiz.columnconfigure(9,weight=1)

raiz.rowconfigure(0,weight=1)
raiz.rowconfigure(1,weight=1)
raiz.rowconfigure(3,weight=1)
raiz.rowconfigure(2,weight=1)
raiz.rowconfigure(4,weight=1)
raiz.rowconfigure(5,weight=1)
raiz.rowconfigure(6,weight=1)

# Cuadros de entrada y salida de texto con scrooll. Etiquetas correspondientes. Excepto porcentajes

Label(frame,text="ADN/ARN",font=("Times New Roman",13)).grid(row=0,column=0,sticky="w",padx=5,pady=10) 

cad_Silvestre=Text(frame,width="11",height="15")
cad_Silvestre.grid(row=1,column=0,sticky=E+W,padx=5,pady=5,rowspan=4,ipadx=4)


scroll_Silv=Scrollbar(frame, command=cad_Silvestre.yview)
scroll_Silv.grid(row=1,column=1,sticky=N+E+S+W,padx=5,pady=5,rowspan=4)
cad_Silvestre.config(yscrollcommand=scroll_Silv.set)


Label(frame,text="Proteína",font=("Times New Roman",13)).grid(row=6,column=0,sticky="w",padx=5,pady=10)
proteina=Text(frame,width="70",height="15")
proteina.grid(row=7,column=0,sticky="w",padx=5,pady=5,rowspan=4)

scroll_prote=Scrollbar(frame,command=proteina.yview)
scroll_prote.grid(row=7,column=1,sticky=N+S+E+W,padx=5,pady=5)
proteina.config(yscrollcommand=scroll_prote.set)

Label(frame,text="Resultados",font=("Times New Roman",13)).grid(row=6,column=3,sticky="w",padx=5,pady=10)
resultados=Text(frame,width="70",height="15")#,textvariable=resultados_salida
resultados.grid(row=7,column=2,sticky=E,padx=5,pady=5,columnspan=4)#,ipadx=240,ipady=110columnspan=5,

scroll_resul=Scrollbar(frame,command=resultados.yview)
scroll_resul.grid(row=7,column=6,sticky=N+S+E+W,padx=5,pady=5)
resultados.config(yscrollcommand=scroll_resul.set)


#FUNCIONES

def generador_cadenas(maximo=2500):     #Genera cadenas de ADN, de longitud aleatoria
    import random as r
    
    resultados.delete(1.0,END)  #Vacía el contenido de Text
    
    num=r.randint(0,maximo)
    bases=["A","T","C","G"]
    cadena=[]
    for n in range (num):
        b=r.choice(bases)
        cadena.append(b)
    adn="".join(cadena)    
    return resultados.insert(END,adn)

def limpiarCadena(contenedor=cad_Silvestre,abc=["A","G","C","T"]): #Limpia las cadenas para que sean ADN, por defecto
    cad_suc=contenedor.get(1.0,"end-1c")    
    cadena_limpia=""
    cad_suc=cad_suc.upper()
    
    for base in cad_suc:
        if base not in abc:
            continue
        else:
            cadena_limpia=base+cadena_limpia
    cadena_limpia=cadena_limpia[::-1]

    return cadena_limpia

def ImprimirCadenaLimpia(): #Imprime el ADN limpio en pantalla.
    resultados.delete(1.0,END)
    return resultados.insert(END,limpiarCadena())


def contarBase(base,texto):#Da el porcentaje del caracter 'base' en el texto.
    try:
        n_base=0
        base.upper()
        texto.upper()
        n_base=texto.count(base)
        
        if base not in ["A","C","T","G","U","AU","GC","AT"]:
            raise ValueError
        
    except ValueError:          #Si el ADN, no solo contiene bases aceptadas, se limpia la cadena.
        texto=limpiarCadena()
        base.upper()
        texto.upper()
        n_base=texto.count(base)
    return round((n_base/len(texto))*100,2)
    
def transcripcion(cad):    #Devuelve el adn tras pasar por la transcripción.             
    complementaria={'A':'T','T':'A','C':'G','G':'C'}
    ADN_ARN={'A':'U','T':'A','C':'G','G':'C'}
    complement=""
    arn=""
        

    for base in cad:
        base_complement=complementaria[base]
        complement=complement+base_complement       #Se hace la cadena complementaria ADN.
            
    complement=complement[::-1]
        
    for base_2 in complement:       #Paso de ADN a ARN.
        base_ARN=ADN_ARN[base_2]
        arn=arn+base_ARN
            
    arn=arn[::-1]
    
    return arn

def transcripcion_ADN():    #Si hay cualquier error se lanza ValueError
                            #El errro será capturado por el except de las funciones principales.
    
    try:
        cadena=cad_Silvestre.get(1.0,"end-1c")
        arn=transcripcion(cadena)
        return arn
      
    except (ValueError,KeyError):
        raise ValueError

def traduccion(cadena):
    
    proteina=""
    
    código_genético=[("UUU","Phe"),("UUC","Phe"),("UUA","Leu"),("UUG","Leu"),("CUU","Leu"),("CUC","Leu"),("CUA","Leu"),("CUG","Leu"),
                     ("AUU","Iso"),("AUC","Iso"),("AUA","Iso"),("AUG","Met"),("GUU","Val"),("GUC","Val"),("GUA","Val"),("GUG","Val"),
                     ("UCU","Ser"),("UCC","Ser"),("UCA","Ser"),("UCG","Ser"),("CCU","Pro"),("CCC","Pro"),("CCA","Pro"),("CCG","Pro"),
                     ("ACU","Thr"),("ACG","Thr"),("ACA","Thr"),("ACC","Thr"),("GCU","Ala"),("GCC","Ala"),("GCA","Ala"),("GCG","Ala"),
                     ("UAU","Tyr"),("UAC","Tyr"),("UAA","Stop"),("UAG","Stop"),("CAU","His"),("CAC","His"),("CAA","Gln"),("CAG","Gln"),
                     ("AAU","Asn"),("AAC","Asn"),("AAA","Lys"),("AAG","Lys"),("GAU","Asp"),("GAC","Asp"),("GAA","Asp"),("GAG","Asp"),
                     ("UGU","Cys"),("UGC","Cys"),("UGA","Stop"),("UGG","Try"),("CGU","Arg"),("CGC","Arg"),("CGA","Arg"),("CGG","Arg"),
                     ("AGU","Ser"),("AGC","Ser"),("AGA","Ser"),("AGG","Ser"),("GGU","Gly"),("GGC","Gly"),("GGA","Gly"),("GGG","Gly")]
    
    pos_i=cadena.find("AUG")#busca el codon de inicio
    arn=cadena[pos_i:]
    
    for i in range(0,len(arn)+3,3):
        
        codon=arn[i:i+3]
        
        for codon_2 in código_genético:
            
            if codon==codon_2[0] and codon_2[1]!="Stop":#comprobar STOP

                    proteina=proteina+codon_2[1]+"-"
                    
            elif codon==codon_2[0] and codon_2[1]=="Stop":

                break
                return proteina
            else:
                continue
    return proteina

def traduccion_ADN():   #Si traduccion devuelve un Error. Se lanza.
                        #Será capturadao pro el except de funciones principales.

    try:
        cadena=cad_Silvestre.get(1.0,"end-1c")
        prote=traduccion(cadena)
        return prote
      

    except (ValueError,IndexError):
        raise ValueError

def ARN_ADN(arn):
    try:
        ARN_ADN={'A':'T','U':'A','G':'C','C':'G'}
        adn=""
        for base in arn:
            adn=adn+ARN_ADN[base]
        return adn
    except KeyError:
        raise ValueError

def complementaria(cadena_molde):
    try:
        complementaria={'A':'T','T':'A','C':'G','G':'C'}
        compl=""
        for c in cadena_molde:
            compl=compl+complementaria[c]
        return compl
    except KeyError:
        raise ValueError
        

 
def estadisticas_CADENA():  #Se vacían todos los Texts
    adn_control=False
    arn_control=False
    
    
    # Porcentajes de las bases nitrogenadas
    resultados.delete(1.0,END)
    p_GC.set("0.0")      #PORCENTAJES
    p_AT.set("0.0")
    p_AU.set("0.0")
    p_U.set("0.0")
    p_C.set("0.0")
    p_G.set("0.0")
    p_A.set("0.0")
    p_T.set("0.0") 
    
    
    adn=cad_Silvestre.get(1.0,"end-1c")
    if adn=="":
        resultados.insert(END,"La cadena de ADN/ARN es una cadena vacía.")
        return None

    try:

        bases=["A","C","T","G","U","AU","GC","AT"]
        

        if adn.isalpha()==True and adn.find("T")!=-1 and adn.find("U")!=-1:
            
           resultados.insert(END,"Error.No se puede determinar si la cadena es ARN o ADN")
           raise ValueError
           return None
       
        elif adn.isalpha()==True and adn.find("U")!=-1:   #Si la cadena introducida no es un texto lanza ValueError
            arn_control=True
            pass
        elif adn.isalpha()==True and adn.find("T")!=-1:
            adn_control=True
            pass

        else:
            raise ValueError

        
        for base in adn:

            if base not in bases:

                break
                raise ValueError
                
            else:

                continue

            
        p_GC.set(contarBase("GC",adn))      #PORCENTAJES
        p_AT.set(contarBase("AT",adn))
        p_AU.set(contarBase("AU",adn))
        p_U.set(contarBase("U",adn))
        p_C.set(contarBase("C",adn))
        p_G.set(contarBase("G",adn))
        p_A.set(contarBase("A",adn))
        p_T.set(contarBase("T",adn))        #Aquí ya se ha comprobado que no es nada diferente de texto y que dentro de eso, son bases aceptadas
                                            #La cadena limpia es:
        resultados.insert(END,"La cadena original es: "+adn+"\n")

        if adn_control==True:
            resultados.insert(END,"\nLa cadena introducida es ácido desoxirribonucleico, es decir, ADN.\n")
        elif arn_control==True:
            resultados.insert(END,"\nLa cadena introducida es ácido ribonucleico, es decir, ARN.\n")
            
        
        
        if adn_control==True:
            
            adn_transcrito=transcripcion_ADN()  #proceso transcripción. Solo si es ADN
            resultados.insert(END,"\nLa cadena de bases nitrogenadas después de pasar por el proceso de transcripción queda: \n"+adn_transcrito)
            
            proteina=traduccion(transcripcion_ADN())

        elif arn_control==True:
            
            proteina=traduccion_ADN()   #proceso traduccion

            
        resultados.insert(END,"\nLa proteína correspondiente a la cadena de bases nitrogenadas introducida es: \n"+proteina)
        
    except ValueError:
        
        p_GC.set("Error")
        p_AT.set("Error")
        p_AU.set("Error")
        p_U.set("Error")
        p_C.set("Error")
        p_G.set("Error")
        p_A.set("Error")
        p_T.set("Error")
        resultados.insert(END,"\nLa cadena no se ha introducido correctamente.")
        
    return None


def estadisticas_PROTEINA():
    encontrado=False
    fin=False
    código_genético=[("UUU","Phe"),("UUC","Phe"),("UUA","Leu"),("UUG","Leu"),("CUU","Leu"),("CUC","Leu"),("CUA","Leu"),("CUG","Leu"),
                     ("AUU","Iso"),("AUC","Iso"),("AUA","Iso"),("AUG","Met"),("GUU","Val"),("GUC","Val"),("GUA","Val"),("GUG","Val"),
                     ("UCU","Ser"),("UCC","Ser"),("UCA","Ser"),("UCG","Ser"),("CCU","Pro"),("CCC","Pro"),("CCA","Pro"),("CCG","Pro"),
                     ("ACU","Thr"),("ACG","Thr"),("ACA","Thr"),("ACC","Thr"),("GCU","Ala"),("GCC","Ala"),("GCA","Ala"),("GCG","Ala"),
                     ("UAU","Tyr"),("UAC","Tyr"),("UAA","Stop"),("UAG","Stop"),("CAU","His"),("CAC","His"),("CAA","Gln"),("CAG","Gln"),
                     ("AAU","Asn"),("AAC","Asn"),("AAA","Lys"),("AAG","Lys"),("GAU","Asp"),("GAC","Asp"),("GAA","Asp"),("GAG","Asp"),
                     ("UGU","Cys"),("UGC","Cys"),("UGA","Stop"),("UGG","Try"),("CGU","Arg"),("CGC","Arg"),("CGA","Arg"),("CGG","Arg"),
                     ("AGU","Ser"),("AGC","Ser"),("AGA","Ser"),("AGG","Ser"),("GGU","Gly"),("GGC","Gly"),("GGA","Gly"),("GGG","Gly")]
    
    prote=proteina.get(1.0,"end-1c")
    
    
    arn=""
    resultados.delete(1.0,END)
    p_GC.set("0.0")      #PORCENTAJES
    p_AT.set("0.0")
    p_AU.set("0.0")
    p_U.set("0.0")
    p_C.set("0.0")
    p_G.set("0.0")
    p_A.set("0.0")
    p_T.set("0.0") 
    
    try:
        if prote=="":
            resultados.insert(END,"No ha introducido ninguna proteína")
            return None
        else:
            
        
        
            lista_codones=prote.split("-")
            
            for codon in lista_codones:
                for triplete in código_genético:
                    if codon==triplete[1]:
                        encontrado=True
                if encontrado==False:
                    raise ValueError
                    
                else:
                    encontrado=False
                
            
            for codon in lista_codones:
                for triplete in código_genético:
                    if codon=="Stop":
                        fin=True
                        break
                    elif codon==triplete[1]:
                        encontrado=True
                        arn=arn+triplete[0]
    
                        break
                    else:
                        continue
                if fin==True:
                    break
                elif encontrado==False:    #significa que el codon de la proteína no se ha encontrado
                    raise ValueError
                    break
                    return None
                else:
                    encontrado=False     #Para poder volver a empezar la búsqueda del siguiente codón
    
            
            resultados.insert(END,"\nLa cadena de ácido ribonucleico mensajero correspondiente es: \n"+arn+"\n")
            p_GC.set(contarBase("GC",arn))      #PORCENTAJES
            p_AT.set(contarBase("AT",arn))
            p_AU.set(contarBase("AU",arn))
            p_U.set(contarBase("U",arn))
            p_C.set(contarBase("C",arn))
            p_G.set(contarBase("G",arn))
            p_A.set(contarBase("A",arn))
            p_T.set(contarBase("T",arn))
            adn=ARN_ADN(arn)
            resultados.insert(END,"\nEl ácido desoxirribonucleico correspondiente a la cadena de ARNm es: \n"+adn+"\n")
            adn_compl=complementaria(adn)
            resultados.insert(END,"\nLa cadena complementaria al ácido desoxirribonucleico del que provienen la proteína es: \n"+adn_compl+"\n")
    
    except ValueError:
        p_GC.set("Error")
        p_AT.set("Error")
        p_AU.set("Error")
        p_U.set("Error")
        p_C.set("Error")
        p_G.set("Error")
        p_A.set("Error")
        p_T.set("Error")
        resultados.insert(END,"La proteina no se ha introducido correctamente")
        
    
  
frame.pack(fill="both",expand=True)                     # Configuración del frame
frame.config(bg="#69f0b6",width="1250",height="655")
frame.config(bd=15,relief="groove")
                                                        
                                                   

#Etiquetas y cuadros de texto de los porcentajes de bases
Label(frame,text="Porcentaje de bases nitrogenadas",font=("Times New Roman",13)).grid(row=0,column=3,sticky="w",padx=5,pady=5,columnspan=3)


Label(frame,text="%GC",font=("Times New Roman",13)).grid(row=1,column=2,sticky=E,padx=5,pady=0)
GC=Entry(frame,font=("Times New Roman",13),textvariable=p_GC,state="readonly")
GC.grid(row=1,column=3,sticky=E,padx="5",pady="5",ipadx=0,ipady=5)

Label(frame,text="%AT",font=("Times New Roman",13)).grid(row=2,column=2,sticky=E,padx=5,pady=0)
AT=Entry(frame,font=("Times New Roman",13),textvariable=p_AT,state="readonly")
AT.grid(row=2,column=3,sticky=E,padx=5,pady=5,ipadx=0,ipady="5")

Label(frame,text="%AU",font=("Times New Roman",13)).grid(row=3,column=2,sticky=E,padx=5,pady=0)
AU=Entry(frame,font=("Times New Roman",13),textvariable=p_AU,state="readonly")
AU.grid(row=3,column=3,sticky=E,padx=5,pady=5,ipadx=0,ipady="5")

Label(frame,text="%G",font=("Times New Roman",13)).grid(row=1,column=4,sticky=E,padx=5,pady=10)
G=Entry(frame,font=("Times New Roman",13),textvariable=p_G,state="readonly")
G.grid(row=1,column=5,sticky="W",padx=5,pady=5,ipadx=0,ipady="5")

Label(frame,text="%A",font=("Times New Roman",13)).grid(row=2,column=4,sticky=E,padx=5,pady=10)
A=Entry(frame,font=("Times New Roman",13),textvariable=p_A,state="readonly")
A.grid(row=2,column=5,sticky="w",padx=5,pady=5,ipadx=0,ipady="5")

Label(frame,text="%T",font=("Times New Roman",13)).grid(row=4,column=4,sticky=E,padx=5,pady=10)
T=Entry(frame,font=("Times New Roman",13),textvariable=p_T,state="readonly")
T.grid(row=4,column=5,sticky="w",padx=5,pady=5,ipadx=0,ipady="5")

Label(frame,text="%C   ",font=("Times New Roman",13)).grid(row=4,column=2,sticky=E,padx=5,pady=10)
C=Entry(frame,font=("Times New Roman",13),textvariable=p_C,state="readonly")
C.grid(row=4,column=3,sticky="w",padx=5,pady=5,ipadx=0,ipady=5)

Label(frame,text="%U",font=("Times New Roman",13)).grid(row=3,column=4,sticky=E,padx=5,pady=10)
U=Entry(frame,font=("Times New Roman",13),textvariable=p_U,state="readonly")
U.grid(row=3,column=5,sticky="w",padx=5,pady=5,ipadx=0,ipady="5")
#Botones
boton_limpiarCadena=B(frame, text="Limpiar ADN/ARN",font=("Times New Roman",13),command=ImprimirCadenaLimpia)
boton_limpiarCadena.grid(row=5,column=2,columnspan=1,padx=5,pady=5,ipadx=4,ipady=5,sticky=S+N+E+W)

boton_estadisticasCad=B(frame,text="Estadisticas ADN",font=("Times New Roman",13),command=estadisticas_CADENA)
boton_estadisticasCad.grid(row=5,column=3,columnspan=1,padx=5,pady=5,ipadx=5,ipady=5,sticky=S+N+E+W)

boton_estadisticasProte=B(frame,text="Estadisticas Proteina",font=("Times New Roman",13),command=estadisticas_PROTEINA)
boton_estadisticasProte.grid(row=5,column=4,columnspan=1,padx=5,pady=5,ipadx=5,ipady=5,sticky=S+N+E+W)

boton_generador=B(frame,text="Generador de ADN",font=("Times New Roman",13),command=generador_cadenas)
boton_generador.grid(row=5,column=5,padx=5,columnspan=1,pady=5,ipadx=5,ipady=5,sticky=S+N+E+W)

raiz.mainloop()

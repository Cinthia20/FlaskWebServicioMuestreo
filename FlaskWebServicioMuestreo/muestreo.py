"""
This script runs the application using a development server.
It contains the definition of routes and views for the application.
"""
import requests
import psycopg2
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri


from flask import Flask, jsonify, request
from FlaskWebServicioMuestreo import app

#app = Flask(__name__)



# Make the WSGI interface available at the top level so wfastcgi can get it.
#wsgi_app = app.wsgi_app

def obtenerPTenRBT(division, zona, anio):
    if division=='0':
        division=''
    if zona=='0':
        zona=''

    #conexion = psycopg2.connect("dbname=PTenBT user=postgres password=Contrasenia2020")
    conexion = psycopg2.connect(host="localhost", database="PTenBT", user="postgres", password="Contrasenia2020")

    # Creamos el cursor con el objeto conexion
    crsPTTotal = conexion.cursor()
    # Ejecutamos una consulta
    sql = "SELECT anio, clavedivision, clavezona, idbancotransformador, ptredes + ptacometidas + ptmedidores as pttotalrbt FROM   ptenbt.fntbl_obtenerptporbanco('" + division + "', '" + zona + "', '', '', " + str(anio) + "); "
    crsPTTotal.execute(sql)

    # Recorremos los resultados y finalmente se depositan en una matriz
    arrPTTotal = [['anio', 'clavedivision', 'clavezona', 'idbancotransformador', 'pttotalrbt']]
    
    for anio, clavedivision, clavezona, idbancotransformador, pttotalrbt in crsPTTotal.fetchall() : 
        valores = [anio, clavedivision, clavezona, idbancotransformador, pttotalrbt]
        arrPTTotal.append(valores)

    # Cerramos la conexión
    conexion.close()

    return arrPTTotal

def deMatrizaLista(matriz):
    lista = []
    for i in range(0, len(matriz)):
        lista.append(matriz[i])
    
    return lista


def obtenerMuestreoDesdeR(datosPT):
    pandas2ri.activate()
    r = robjects.r
    r.assign('datos', pandas2ri.py2rpy(datosPT))

    r(''' 
        library(survey)

        # PARÁMETRO DE INTERÉS: TOTAL
        Tot = sum(datos$pttotalrbt)

        # PARÁMETROS DE LA SIMULACIÓN
        ns = seq(100,1500,50)
        N = nrow(datos)
        M = 500

        # 1) MUESTREO ALEATORIO SIMPLE
        set.seed(9654)
        total <- 0
        totales <- matrix(0,nrow=M,ncol=length(ns))
        for (i in 1:length(ns)) {
            for(j in 1:M) { 
                muestra = sample(datos$pttotalrbt,ns[i])
                total = N*mean(muestra)
                totales[j,i] = total
            }
        }

        IC <- matrix(0,nrow=length(ns),ncol=2)
        ancho <- 0
        for(i in 1:length(ns)) {
            IC[i,1] = quantile(totales[,i],0.025)
            IC[i,2] = quantile(totales[,i],0.975)
        }

        # 2) MUESTREO ESTRATIFICACIÓN
        # CONSTRUCCIÓN DE LOS ESTRATOS
        datos02 <- datos[order(datos$pttotalrbt),]

        perdidas = datos02$pttotalrbt/Tot
        cumperd = cumsum(perdidas)
        datos02 <- data.frame(datos02, perdAcum = cumperd)

        estr01 <- datos02[which(datos02$perdAcum <= 0.20),]
        estr02 <- datos02[which(datos02$perdAcum > 0.20 & datos02$perdAcum <= 0.40),]
        estr03 <- datos02[which(datos02$perdAcum > 0.40 & datos02$perdAcum <= 0.80),]
        estr04 <- datos02[which(datos02$perdAcum > 0.80),]

        e01 = data.frame(estr01,estrato=rep("e01",nrow(estr01)))
        e02 = data.frame(estr02,estrato=rep("e02",nrow(estr02)))
        e03 = data.frame(estr03,estrato=rep("e03",nrow(estr03)))
        e04 = data.frame(estr04,estrato=rep("e04",nrow(estr04)))

        datos02 <- rbind(e01,e02,e03,e04)

        # TAMAÑO DE LOS ESTRATOS
        Nh = table(datos02$estrato)

        # NÚMERO DE ESTRATOS
        L = 4

        # TAMAÑO DE ESTRATOS ACUMULADOS
        Nhacumulado = 0
        suma = 0
        for(i in 1:L)
        {
            suma = suma + Nh[i]
            Nhacumulado[i] = suma 
        }

        # INICIALIZAMOS SEMILLA
        set.seed(5580)

        # IC2 ES LA MATRIZ DONDE GUARDAMOS LOS INTERVALOS DE CONFIANZA
        IC2 <- matrix(0,nrow=length(ns),ncol=2)

        for(l in 1:length(ns)) {
            n <- ns[l]
            # ASIGNACIÓN DE NEYMAN
            s = as.numeric(sqrt(by(datos02$pttotalrbt,datos02$estrato,var)))
            nh <- round(n*Nh*s/sum(Nh*s),0)

            # FACTORES DE CORRECCIÓN POR POBLACIÓN FINITA
            fpcdis <- matrix(0,nrow=sum(nh),ncol=1)
            k = 1
            for(i in 1:L) {
                for(j in 1:nh[i]) {
                    fpcdis[k] <- Nh[i]
                    k <- k + 1
                } 
            }

            resultado <- matrix(0,nrow=M,ncol=1)
            for (i in 1:M) {
                # VECTOR QUE INDICA LOS INDICES SELECCIONADOS
                hileras <- c(
                    sample(Nhacumulado[1],nh[1]),
                    sample((Nhacumulado[1]+1):Nhacumulado[2],nh[2]),
                    sample((Nhacumulado[2]+1):Nhacumulado[3],nh[3]),
                    sample((Nhacumulado[3]+1):Nhacumulado[4],nh[4])
                )
                muestra <- datos02[hileras,]
                muestra_estratificada = data.frame(muestra,fpcdis)

                # CONSTRUIMOS OBJETO DE CLASE MUESTRA ESTRATIFICADA
                muestra.str <- svydesign(id=~1,strata=~estrato,data=muestra_estratificada,fpc=fpcdis)

                # ESTIMAMOS EL TOTAL
                resultado[i] = svytotal(~pttotalrbt,muestra.str, deff=TRUE)[1]
            }

            IC2[l,1] = quantile(resultado,0.025)
            IC2[l,2] = quantile(resultado,0.975)
        }   
        
    ''')

    ptTotal = r('Tot')
    print(ptTotal)
    simpleInferior = r('c(IC[,1])')
    simpleSuperior = r('c(IC[,2])')
    estratificadoInferior = r('c(IC2[,1])')
    estratificadoSuperior = r('c(IC2[,2])')
    ejeX = r('c(ns)')

    lista = []
    for i in range(0, len(ejeX)):
        renglon = [ejeX[i], ptTotal[0], simpleInferior[i], simpleSuperior[i], estratificadoInferior[i], estratificadoSuperior[i]]
        lista.append(renglon)
    

    return jsonify({
        "Descripcion": ["col0: Eje X", 
                        "col1: Perdida Tecnica Total en RBT"
                        "col2: Muestreo Aleatorio Simple Inferior",
                        "col3: Muestreo Aleatorio Simple Superior", 
                        "col4: Muestreo Estratificado Inferior",
                        "col5: Muestreo Estratificado Superior"],
        "datosGrafica": lista,
        })


def obtenerMuestreo(division, zona, anio):
  # ---------------------------- Main ---------------------

  #Se obtienen los datos en forma de matriz
  dtPTenRBT = obtenerPTenRBT(division, zona, anio)

  # Se genera el dataframe
  dfPT = pd.DataFrame(data=dtPTenRBT, columns=dtPTenRBT.pop(0))

  #se ejecuta desde R y en la misma función se hace el pblado de las variables globales
  #return obtenerMuestreoDesdeR(dfPT)

  marco = obtenerMuestreoDesdeR(dfPT)

  return marco
  


@app.route('/')
def hello():
    return  "Para obtner datos proporcione anio, [division], [zona]"+ 'http://localhost:60126/muestreo/DB/01/2019'

# agregar un nuevo producto
#@app.route('/muestreo/', methods=['POST'])
#def addProduct():
#    division = request['DB']   
#    zona = request['zona']
#    anio = request['anio']   

#    return obtenerMuestreo(division, zona, anio)

@app.route('/muestreo/<division>/<zona>/<anio>')
#@app.route('/muestreo/division='+'<division>'+'&zona=<zona>'+'&anio=<anio>')
def addProduct(division,zona,anio):
    #msg="Hello "+division+","+zona+","+anio+"!!! Welcome to my Website"
    #return msg     
    return obtenerMuestreo(division,zona, anio)


##if __name__ == '__main__':
# #   import os
#    HOST = os.environ.get('SERVER_HOST', 'localhost')
#    try:
#        PORT = int(os.environ.get('SERVER_PORT', '5555'))
#    except ValueError:
#        PORT = 5555
#    app.run(HOST, PORT, debug =True)
#    #app.run(HOST, PORT, debug=False)
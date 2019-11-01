using LinearAlgebra;
using Plots;

#Función que nos regresará la matriz "Skew Symmetric" dado un arreglo de 3 entradas.
function skew_Symm(A::Array)
	return reshape([0, A[3], -A[2], -A[3], 0, A[1], A[2], -A[1], 0],3,3);
end

#Función que dado un vector de posición de la partícula y la longitud de las aristas de la caja cúbica (centrada en el origen) que contiene a la partícula, nos regresa la posición de la 
#partícula dentro de la caja tras aplicarle condiciones periódicas a la frontera.
#Posicion: Arreglo con las coordenadas de la partícula.
#LadoCaja: Longitud de la arista de la caja.
function caja_Periodica!(Posicion::Array, LadoCaja::Float64)
	N = length(Posicion); #Dimensión del espacio donde vive la partícula
	Caja = floor.((Posicion .+ LadoCaja/2.) ./ LadoCaja) * LadoCaja; #Nos indica la coordenada de la caja en la cual se encuentra la partícula antes de las condiciones periódicas

	for i in 1:N #Aplicamos condiciones periódicas a la frontera a cada una de las coordenadas de la partícula
		if (-LadoCaja/2. > Posicion[i]) || (Posicion[i] >= LadoCaja/2.) #Condición de que la partícula se ha salido de la caja 
			Posicion[i] = mod(Posicion[i] + LadoCaja/2., LadoCaja) - LadoCaja/2.; #Obtenemos la coordenada relativa al centro de la caja que contiene a la partícula
		end
	end

	return Posicion, Caja
end

#Función que regresa la matriz de rotación que transforma al vector (0,0,1) en el vector dado por el usuario.
#Vector: Arreglo con las coordenadas del vector dado por el usuario.
function mat_Rot(Vector::Array)
	E1 = [0.0, 0.0, 1.0]; #Vector unitario en Z
	Vector = Vector/norm(Vector); #Normalicemos el "Vector" dado por el usuario
	Coseno = dot(E1,Vector); #Producto punto entre los vectores E1 y Vector
    (Coseno == -1) ? throw(ArgumentError("Brinde un vector que no sea colineal al vector (0,0,1)")) : nothing; #Verificamos que "Vector" no sea colineal al vector (0,0,1)
    Cruz = cross(E1,Vector); #Vector generado por el producto cruz de los vectores E1 y Vector
	Discriminante = 1/(1 + Coseno);

	return Matrix{Float64}(I,3,3) + skew_Symm(Cruz) + Discriminante*(skew_Symm(Cruz)*skew_Symm(Cruz)); #Obtenemos la matriz de rotación que lleva (0,0,1) -> Vector
end

#Función que genera un vector unitario aleatorio con una distribución equiproble sobre el disco unitario.
#Mat_Z_Cil: Arreglo de 3x3 correspondiente a la matriz que lleva (0,0,1) -> Vector, con Vector un vector arbitrario en el espacio.
function vel_Rand(Mat_Z_Cil::Array)
	θ = (2*π)*rand(); #Generamos un ángulo entre 0 y 2*Pi arbitrario
	Vel_Sin_Girar = [cos(θ), sin(θ), 0]; #Generamos un vector sobre el plano XY que forma un ángulo respecto al eje X igual al ángulo generado en la línea previa.
	return Mat_Z_Cil*Vel_Sin_Girar; #Rotamos el vector velocidad en el plano XY para que sea ahora ortogonal al eje del vector "Vector" que define a "Mat_Z_Cil"
end

#Función que regresa los centros de la caja principal [centrada en (0,0,0)] y de las 26 cajas adyacentes a ella.
#LadoCaja: Longitud de la arista de la caja.
function centros_Cajas_3D(LadoCaja::Float64)
	A = [[-LadoCaja, LadoCaja, -LadoCaja], 
        [0.0, LadoCaja, -LadoCaja], 
        [LadoCaja, LadoCaja, -LadoCaja], 
        [-LadoCaja, 0.0, -LadoCaja], 
        [0.0, 0.0, -LadoCaja], 
        [LadoCaja, 0.0, -LadoCaja], 
        [-LadoCaja, -LadoCaja, -LadoCaja], 
        [0.0, -LadoCaja, -LadoCaja], 
        [LadoCaja, -LadoCaja, -LadoCaja], 
        [-LadoCaja, LadoCaja, 0.0], 
        [0.0, LadoCaja, 0.0], 
        [LadoCaja, LadoCaja, 0.0], 
        [-LadoCaja, 0.0, 0.0], 
        [0.0, 0.0, 0.0], 
        [LadoCaja, 0.0, 0.0], 
        [-LadoCaja, -LadoCaja, 0.0], 
        [0.0, -LadoCaja, 0.0], 
        [LadoCaja, -LadoCaja, 0.0], 
        [-LadoCaja, LadoCaja, LadoCaja], 
        [0.0, LadoCaja, LadoCaja], 
        [LadoCaja, LadoCaja, LadoCaja], 
        [-LadoCaja, 0.0, LadoCaja], 
        [0.0, 0.0, LadoCaja], 
        [LadoCaja, 0.0, LadoCaja], 
        [-LadoCaja, -LadoCaja, LadoCaja], 
        [0.0, -LadoCaja, LadoCaja], 
        [LadoCaja, -LadoCaja, LadoCaja]]
    return A
end

#Función que calcula la longitud que debe tener el cilindro que insertaremos dentro de la caja cúbica central.
#LadoCaja: Longitud de la arista de la caja.
#Mat_Cil_Z: Matriz que transforma "Vector" -> (0,0,1) para algún "Vector" dado por el usuario.
function longitud_Cilindro(LadoCaja::Float64, Mat_Cil_Z::Array)
	#Definimos las coordenadas de las esquinas del cubo sin rotar.
    Arreglo_Aristas_No_Rot  =   [[-(LadoCaja)/2, (LadoCaja)/2, -(LadoCaja)/2], 
                                [-(LadoCaja)/2, -(LadoCaja)/2, -(LadoCaja)/2], 
                                [(LadoCaja)/2, (LadoCaja)/2, -(LadoCaja)/2], 
                                [(LadoCaja)/2, -(LadoCaja)/2, -(LadoCaja)/2],
                                [-(LadoCaja)/2, (LadoCaja)/2, (LadoCaja)/2], 
                                [-(LadoCaja)/2, -(LadoCaja)/2, (LadoCaja)/2], 
                                [(LadoCaja)/2, (LadoCaja)/2, (LadoCaja)/2], 
                                [(LadoCaja)/2, -(LadoCaja)/2, (LadoCaja)/2]];

    #Giramos el cubo de modo que el cilindro quede alineado con el eje Z, esto nos da las siguientes coordenadas.
    Arreglo_Aristas_Rot = []; #Arreglo donde colocaremos las aristas del cubo tras rotar el espacio "Vector" -> (0,0,1)

    for i in Arreglo_Aristas_No_Rot #Tomamos las coordenadas de cada una de las aristas del cubo sin rotar
    	x = Mat_Cil_Z*i; #Rotamos el vector asociado a la arista "i"
    	push!(Arreglo_Aristas_Rot, x);
    end

    #Obtenemos la máxima longitud entre cada uno de los puntos en su coordenada en Z
    Arreglo_Longitud_Z = Float64[];

    for i in 1:(length(Arreglo_Aristas_Rot)-1)
    	for j in (i+1):length(Arreglo_Aristas_Rot)
    		LongitudZ = abs(Arreglo_Aristas_Rot[i][3] - Arreglo_Aristas_Rot[j][3]); #Dados los vértices i y j del cubo rotado, tomamos la diferencia entre las alturas de ambos.
    		push!(Arreglo_Longitud_Z, LongitudZ);
    	end
    end

    return maximum(Arreglo_Longitud_Z); #Obtenemos la máxima separación entre la altura de las aristas del cubo rotado.
end

#Función que regresa un arreglo con un determinado número de puntos que estén dentro de la caja central y fuera del cilindro (con sus subcilindros por condiciones periódicas a la frontera).
#RadioCil: Flotante con el radio del cilindro que insertamos en la caja cúbica central.
#LadoCil: Longitud del cilindro
#LadoCaja: Longitud de la arista de la caja.
#Mat_Cil_Z: Matriz de rotación EjeCil -> (0,0,1).
#Iteraciones: Número de puntos máximo que se desea obtener tales que estén dentro de la caja central y fuera del cilindro (con sus subcilindros por condiciones periódicas a la frontera).
function posiciones_Iniciales(RadioCil::Float64, LongitudCil::Float64, LadoCaja::Float64, Mat_Cil_Z::Array, Iteraciones::Int64)
	Arreglo_Centros_Cilindros = centros_Cajas_3D(LadoCaja); #Definimos los centros de todas las cajas adyacentes a la caja central.
	Arreglo_Centros_Cilindros_Rotados = []; #Arreglo que contendrá los centros de las 27 cajas tras rotar el espacio para que EjeCil -> (0,0,1)

    #Obtengamos estos centros tras rotar el espacio "EjeCil -> (0,0,1)" y tener los cilindros verticales
	for i in Arreglo_Centros_Cilindros
		x = Mat_Cil_Z*i; #Aplicamos a cada uno de los vértices la matriz de rotación "EjeCil -> (0,0,1)"
		push!(Arreglo_Centros_Cilindros_Rotados, x);
	end

	Arreglo_Posiciones = []; #Arreglo donde irán las coordenadas de puntos válidos

	for i in 1:Iteraciones
		#Llaves para iniciar el código
		Calcular = true;
		Intentos = 0;

		while Calcular == true && Intentos < 1e4
			Intentos += 1;

			#En este punto coloquemos una partícula, sin importarnos si estamos dentro o fuera del cilindro
			PrevX = rand() - 0.5;
			PrevY = rand() - 0.5;
			PrevZ = rand() - 0.5;
			Prev_Pos = [PrevX, PrevY, PrevZ]; #Posición de la partícula dentro de la caja central (sin rotar)

			#Veamos si ese punto es un punto válido.
            Pos = Mat_Cil_Z*Prev_Pos; #Obtenemos las coordenadas de la partícula tras rotar el espacio "EjeCil -> (0,0,1)"

            Dentro_Cilindro = false; #Indicador de si la partícula cayó dentro del cilindro o no

			#Ahora coloquemos el origen en el centro de cada cilindro rotado y veamos si el punto cae dentro.
			for i in Arreglo_Centros_Cilindros_Rotados
				Posicion_Punto = Pos - i; #Punto de la partícula a prueba tras la rotación y traslación espacial del cilindro i-ésimo al origen.
				(-(LongitudCil/2) < Posicion_Punto[3] < LongitudCil/2) && (((Posicion_Punto[1])^2 + (Posicion_Punto[2])^2) < (RadioCil)^2) ? (Dentro_Cilindro = true; break) : (nothing);
            end
                
            #Si el punto propuesto resulta que cae dentro de un cilindro, lo descartamos, si no, lo aceptamos
            Dentro_Cilindro ? (nothing) : (push!(Arreglo_Posiciones, Prev_Pos); Calcular = false);
		end
	end

	return Arreglo_Posiciones
end

#Función que regresa un arreglo con un determinado número de puntos que estén dentro de la caja central y dentro del cilindro (con sus subcilindros por condiciones periódicas a la frontera).
#RadioCil: Flotante con el radio del cilindro que insertamos en la caja cúbica central.
#LadoCil: Longitud del cilindro
#LadoCaja: Longitud de la arista de la caja.
#Mat_Cil_Z: Matriz de rotación EjeCil -> (0,0,1).
#Iteraciones: Número de puntos máximo que se desea obtener tales que estén dentro de la caja central y dentro del cilindro (con sus subcilindros por condiciones periódicas a la frontera).
function posiciones_Iniciales_2(RadioCil::Float64, LongitudCil::Float64, LadoCaja::Float64, Mat_Cil_Z::Array, Iteraciones::Int64)
	Arreglo_Centros_Cilindros = centros_Cajas_3D(LadoCaja); #Definimos los centros de todas las cajas adyacentes a la caja central.
	Arreglo_Centros_Cilindros_Rotados = []; #Arreglo que contendrá los centros de las 27 cajas tras rotar el espacio para que EjeCil -> (0,0,1)

    #Obtengamos estos centros tras rotar el espacio "EjeCil -> (0,0,1)" y tener los cilindros verticales
	for i in Arreglo_Centros_Cilindros
		x = Mat_Cil_Z*i; #Aplicamos a cada uno de los vértices la matriz de rotación "EjeCil -> (0,0,1)"
		push!(Arreglo_Centros_Cilindros_Rotados, x);
	end

	Arreglo_Posiciones = []; #Arreglo donde irán las coordenadas de puntos válidos

	for i in 1:Iteraciones
		#Llaves para iniciar el código
		Calcular = true;
		Intentos = 0;

		while Calcular == true && Intentos < 1e4
			Intentos += 1;

			#En este punto coloquemos una partícula, sin importarnos si estamos dentro o fuera del cilindro
			PrevX = rand() - 0.5;
			PrevY = rand() - 0.5;
			PrevZ = rand() - 0.5;
			Prev_Pos = [PrevX, PrevY, PrevZ]; #Posición de la partícula dentro de la caja central (sin rotar)

			#Veamos si ese punto es un punto válido.
            Pos = Mat_Cil_Z*Prev_Pos; #Obtenemos las coordenadas de la partícula tras rotar el espacio "EjeCil -> (0,0,1)"

            Dentro_Cilindro = false; #Indicador de si la partícula cayó dentro del cilindro o no

			#Ahora coloquemos el origen en el centro de cada cilindro rotado y veamos si el punto cae dentro.
			for i in Arreglo_Centros_Cilindros_Rotados
				Posicion_Punto = Pos - i; #Punto de la partícula a prueba tras la rotación y traslación espacial del cilindro i-ésimo al origen.
				(-(LongitudCil/2) < Posicion_Punto[3] < LongitudCil/2) && (((Posicion_Punto[1])^2 + (Posicion_Punto[2])^2) < (RadioCil)^2) ? (Dentro_Cilindro = true; break) : (nothing);
            end
                
            #Si el punto propuesto resulta que cae dentro de un cilindro, lo descartamos, si no, lo aceptamos
            Dentro_Cilindro ? (push!(Arreglo_Posiciones, Prev_Pos); Calcular = false) : (nothing);
		end
	end

	return Arreglo_Posiciones
end

#Función que grafica los cilindros que quedan dentro de la caja central.
#Mat_Z_Cil: Matriz de rotación (0,0,1) -> EjeCil.
#RadioCil: Radio del cilindro que introducimos a la caja principal.
#LongitudCil: Longitud del cilindro que introducimos a la caja principal.
#LadoCaja: Longitud de una de las aristas de la caja cúbica centrada en el origen.
#Caja_i_X: Coordenada X del centro de la caja i-ésima donde está contenido el cilindro que queremos dibujar
#Caja_i_Y: Coordenada Y del centro de la caja i-ésima donde está contenido el cilindro que queremos dibujar
#Caja_i_Z: Coordenada Z del centro de la caja i-ésima donde está contenido el cilindro que queremos dibujar
function dibuja_Cilindros(Mat_Z_Cil::Array, RadioCil::Float64, LongitudCil::Float64, LadoCaja::Float64, Caja_i_X::Float64, Caja_i_Y::Float64, Caja_i_Z::Float64)
	ZZ = range(-LongitudCil/2., stop=LongitudCil/2., length=100); #Definimos el número de capas que tomaremos a lo largo de la altura del cilindro.
	θ = range(0.0, stop=2*π, length=10); #Definimos el número de ángulos a considerar en cada una de las capas que tomaremos.

	X = Float64[]; #Arreglo donde irán las coordenadas X de los puntos del cilindro que estamos considerando que caen dentro de la caja central.
	Y = Float64[]; #Arreglo donde irán las coordenadas Y de los puntos del cilindro que estamos considerando que caen dentro de la caja central.
	Z = Float64[]; #Arreglo donde irán las coordenadas Z de los puntos del cilindro que estamos considerando que caen dentro de la caja central.

	for j in θ
		for i in ZZ
			jj = (2*π)*rand(); #Generamos un ángulo arbitrario entre 0 y 2*Pi
			ii = rand(ZZ); #Tomamos una altura arbitraria entre -LongitudCil/2 y LongitudCil/2
			xx = [RadioCil*cos(jj); RadioCil*sin(jj); ii]; #Generamos un punto del cilindro (vertical) que esté a la altura considerada y en el ángulo considerado
			xx = Mat_Z_Cil*xx; #Le aplicamos la rotación que lleva "(0,0,1) -> EjeCil" al punto previo
			xx = xx + [Caja_i_X, Caja_i_Y, Caja_i_Z]; #Le sumamos al punto previo las coordenadas de la caja que contiene al cilindro que estamos considerando

			if -LadoCaja/2. < xx[1] < LadoCaja/2. && -LadoCaja/2. < xx[2] < LadoCaja/2. && -LadoCaja/2. < xx[3] < LadoCaja/2. #Si el punto cae dentro de la caja central lo guardamos
				push!(X, xx[1]);
				push!(Y, xx[2]);
				push!(Z, xx[3]);
			end
		end
	end

	#Graficamos la superficie 3D generada por los puntos que quedaron dentro de la caja central
	A = scatter(X,Y,Z, aspect_ratio = :equal, markersize = 1, key = false);

	return A
end

#Función que calcula los dos puntos de intersección entre dos círculos, uno centrado en (0,0) y otro centrado en (X,0).
#Radio1: Radio del círculo centrado en (0,0).
#Radio2: Radio del círculo centrado en (X,0).
#X: Coordenada X del segundo círculo.
function interseccion_1(Radio1::Float64, Radio2::Float64, X::Float64)	
	X1 = (Radio1^2 - Radio2^2 + X^2)/(2*X); #Coordenada en X de la intersección de los dos círculos (si es que existe)
	Discriminante = -Radio1^4 + 2*(Radio1^2)*(Radio2^2) - Radio2^4 + 2*(Radio1^2)*(X^2) + 2*(Radio2^2)*(X^2) - X^4; #Valor que nos indica si los círculos se intersectan o no

	if Discriminante >= 0.0 #En caso de ser cierta la desigualdad, hay intersección de los círculos
		Y1 = -sqrt(Discriminante)/(2*X); #Coordenada en Y de la intersección de los dos círculos
		return [X1; Y1], [X1; -Y1] #Regresamos las coordenadas de los dos puntos de intersección de los círculos
	else
		return [Inf; Inf], [Inf; Inf] #Caso contrario, no hay intersección y regresamos como puntos de intersección el infinito
	end
end

#Función que calcula la matriz de rotación V1 -> V2, siendo V1 un vector y V2 otro vector en un espacio 2D. De igual forma regresa la "Distancia" entre los puntos.
#V1: Arreglo con las coordenadas (X,Y) del vector V1.
#V2: Arreglo con las coordenadas (X,Y) del vector V2.
function rotacion(V1::Array, V2::Array)
	Distancia = norm(V1-V2); #Distancia entre los puntos V1 y V2

	CC = dot(V2 - V1, [1.0; 0.0])/Distancia; #Cos(Theta), siendo Theta el ángulo entre el eje X y el vector V2-V1
	SS = dot(V2 - V1, [0.0; 1.0])/Distancia; #Sin(Theta), siendo Theta el ángulo entre el eje X y el vector V2-V1

	return [[CC -SS]; [SS CC]], Distancia
end

#Función que calcula los dos puntos de intersección entre dos círculos, uno con centro en Centro1 y Radio1, otro con centro en Centro2 y Radio2.
#Centro1: Arreglo con las coordenadas (X1,Y1) del centro del círculo 1.
#Radio1: Radio del círculo centrado en (X1,Y1).
#Centro2: Arreglo con las coordenadas (X2,Y2) del centro del círculo 2.
#Radio2: Radio del círculo centrado en (X2,Y2).
function interseccion_2(Radio1::Float64, Radio2::Float64, Centro1::Array, Centro2::Array)
	#Solucionamos la intersección como si los círculos estuvieran en (0,0) y en (d,0).
	Matriz_Rotacion, Distancia = rotacion(Centro1, Centro2);
	Prev_Colision1, Prev_Colision2 = interseccion_1(Radio1, Radio2, Distancia);

	if Prev_Colision1[1] == Inf #Si la colisión se da en infinito, eso no se modificará tras rotar y trasladar
		return [Inf; Inf], [Inf; Inf]
	end

	#Con la matriz de rotación regresamos las soluciones al caso en el que no están centrados los círculos en el origen y en el (d,0).
	Colision1 = Matriz_Rotacion*Prev_Colision1 + Centro1;
	Colision2 = Matriz_Rotacion*Prev_Colision2 + Centro1;

	return Colision1, Colision2
end

#Función que nos regresa el tiempo que tarda en desplazarse la partícula desde "Posicion_Inicial" hasta "Posicion_Final" siguiendo una trayectoria circular con centro en "Centro" y "Radio".
#El tiempo corresponde al ángulo que hay entre "Posicion_Final" y "Posicion_Inicial" considerando que la circunferencia es recorrida en el sentido de las manecillas del reloj.
#Posicion_Inicial: Arreglo con las coordenadas (Xi,Yi) de la posición inicial de la partícula.
#Posicion_Final: Arreglo con las coordenadas (Xf,Yf) de la posición final de la partícula.
#Centro: Arreglo con las coordenadas (X,Y) del centro de la circunferencia que describe la partícula.
#Radio: Radio de la circunferencia con centro en (X,Y) que describe la partícula.
function tiempo_Trayectoria(Posicion_Inicial::Array, Posicion_Final::Array, Centro::Array, Radio::Float64)
	Coseno = (Posicion_Inicial[1] - Centro[1])/Radio; #Coseno del ángulo entre el eje X y la posición inicial de la partícula relativa al centro de su circunferencia
	Seno = (Posicion_Inicial[2] - Centro[2])/Radio; #Seno del ángulo entre el eje X y la posición inicial de la partícula relativa al centro de su circunferencia

	Cosωt = Coseno*((Posicion_Final[1] - Centro[1])/Radio) + Seno*((Posicion_Final[2] - Centro[2])/Radio); #Cos(Theta_f - Theta_i)
	Senωt = Coseno*((Posicion_Final[2] - Centro[2])/Radio) - Seno*((Posicion_Final[1] - Centro[1])/Radio); #Sin(Theta_f - Theta_i)

	if Cosωt > 1.0
		Cosωt = 1.0;
		return 0.0
	end

	if Cosωt < -1.0
		Cosωt = -1.0;
	end

	if Senωt < 0.0
		Tiempo = acos(Cosωt)*Radio; #El tiempo de vuelo corresponde a (Theta_f - Theta_i)*Radio
		return Tiempo
	end

	if Senωt > 0.0
		Tiempo = (2.0*π - acos(Cosωt))*Radio; #El tiempo de vuelo corresponde a [2*Pi - (Theta_f - Theta_i)]*Radio
		return Tiempo
	end
end

#Función que nos regresa el centro de la circunferencia de radio "Radio", que pasa por "Posicion" y con velocidad lineal en "Posicion" igual a "Velocidad".
#Velocidad: Arreglo con coordenadas (Vx,Vy) del vector tangente a la circunferencia paralelo a la velocidad lineal de la partícula en el punto "Posicion".
#Posicion: Arreglo con coordenadas (X,Y) de la posición por la cual pasa la circunferencia deseada.
#Radio: Radio de la circunferencia que describe la partícula.
function centro_Circunferencia(Velocidad::Array, Posicion::Array, Radio::Float64)
	Velocidad = Velocidad/norm(Velocidad); #Vector unitario tangente a la circunferencia con dirección paralela a la velocidad lineal de la partícula
	Velocidad_Ortogonal = [Velocidad[2], -Velocidad[1]]; #Vector unitario que apunta, desde la posición de la partícula "Posicion", al centro de la circunferencia
	Centro = Velocidad_Ortogonal*Radio + Posicion; #Arreglo con las coordenadas del centro de la circunferencia
	return Centro
end

#Funcion que nos da un arreglo con las coordenadas (Vx,Vy) de la velocidad lineal de la partícula en el punto "Posicion" considerando que sigue una trayectoria circular en el sentido de las manecillas del reloj
#centrada en "Centro".
#Posicion: Arreglo con las coordenadas (X,Y) asociadas a la posición de la partícula en una trayectoria circular.
#Centro: Arreglo con las coordenadas (Cx, Cy) asociadas al centro de la circunferencia que describe la partícula al moverse.
function velocidad_Particula(Posicion::Array, Centro::Array)
	Direccion = Posicion - Centro; #Vector que apunta del centro a la posición de la partícula
	Direccion = Direccion/norm(Direccion); #Volvemos unitario el vector previo
	return [Direccion[2], -Direccion[1]] #Regresamos al usuario el vector unitario que apunta en la dirección de la velocidad lineal de la partícula
end

#Función que determina la "Posicion_Final" y la "Velocidad_Final" de una partícula cuya "Posicion_Inicial" y "Velocidad_Inicial" son conocidas. La partícula se desplaza describiendo una circunferencia de radio 
#"Radio" en el sentido de las manecillas del reloj. El tiempo de vuelo es Δt.
#Posicion_Inicial: Arreglo con las coordenadas (Xi, Yi) asociadas a la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las coordenadas (Vxi, Vyi) asociadas a la velocidad inicial de la partícula.
#Radio: Radio de la circunferencia descrita por el desplazamiento de la partícula.
#Δt: Intervalo de tiempo en el que se desea que se desplace la partícula.
function avance(Posicion_Inicial::Array, Velocidad_Inicial::Array, Radio::Float64, Δt::Float64)
	Centro = centro_Circunferencia(Velocidad_Inicial, Posicion_Inicial, Radio); #Localizamos el centro de la circunferencia que describe la partícula al desplazarse
	ω = norm(Velocidad_Inicial)/Radio; #Determinamos la velocidad angular de la partícula
	Posicion_Final = zeros(Float64, 2); #Definimos el arreglo que contendrá las coordenadas de la "Posicion_Final" de la partícula

	Coseno = (Posicion_Inicial[1] - Centro[1])/Radio; #Coseno del ángulo entre el eje X y la posición inicial de la partícula relativa al centro de su circunferencia
	Seno = (Posicion_Inicial[2] - Centro[2])/Radio; #Seno del ángulo entre el eje X y la posición inicial de la partícula relativa al centro de su circunferencia

	Posicion_Final[1] = Radio*(cos(-ω*Δt)*Coseno - sin(-ω*Δt)*Seno) + Centro[1]; #Radio*Cos(Theta_i - w*Δt) + Cx
	Posicion_Final[2] = Radio*(cos(-ω*Δt)*Seno + sin(-ω*Δt)*Coseno) + Centro[2]; #Radio*Sin(Theta_i - w*Δt) + Cy

	Velocidad_Final = velocidad_Particula(Posicion_Final, Centro); #Calculamos la "Velocidad_Final" lineal de la partícula

	return Posicion_Final, Velocidad_Final
end

#Función que calcula la posición de colisión entre la partícula y un obstáculo, así como el tiempo que tarda en realizarse dicha colisión.
#Posicion_Inicial: Arreglo con las coordenadas (Xi, Yi) asociadas a la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las coordenadas (Vxi, Vyi) asociadas a la velocidad inicial de la partícula.
#Radio_Particula: Radio de la circunferencia descrita por el desplazamiento de la partícula.
#Centro_Obstaculo: Arreglo con las coordenadas (Cx, Cy) asociadas al centro del obstáculo con el que la partícula podría colisionar.
#Radio_Obstaculo: Radio del obstáculo circular con el que la partícula podría colisionar
function colision(Posicion_Inicial::Array, Velocidad_Inicial::Array, Radio_Particula::Float64, Centro_Obstaculo::Array, Radio_Obstaculo::Float64)
	Centro_Particula = centro_Circunferencia(Velocidad_Inicial, Posicion_Inicial, Radio_Particula); #Calculamos el centro de la circunferencia que describe la partícula al moverse
	Colision1, Colision2 = interseccion_2(Radio_Particula, Radio_Obstaculo, Centro_Particula, Centro_Obstaculo); #Obtenemos las coordenadas de las posibles colisiones entre la partícula y el obstáculo

	if Colision2[1] == Inf || Colision2[1] == -Inf || Colision2[1] == NaN
		return [Inf; Inf], Inf #Si la colisión no ocurre, regresa al usuario que la colisión se da en infinito en un tiempo infinito
	end

	Tiempo = tiempo_Trayectoria(Posicion_Inicial, Colision2, Centro_Particula, Radio_Particula); #Obtenemos el tiempo que tarda en producirse la colisión.
	return Colision2, Tiempo
end

#Función que calcula la velocidad lineal de la partícula tras colisionar con algún obstáculo.
#Centro_Obstaculo: Arreglo con las coordenadas (Cx, Cy) asociadas al centro del obstáculo con el que la partícula colisionó.
#Posicion_Colision: Arreglo con las coordenadas (X,Y) asociadas a la posición de la partícula al colisionar con el obstáculo.
#Velocidad_Colision: Arreglo con las coordenadas (Vx, Vy) asociadas a la velocidad de la partícula al colisionar con el objetivo.
function velocidad_Tras_Colision(Centro_Obstaculo::Array, Posicion_Colision::Array, Velocidad_Colision::Array)
	Normal_Colision = (Centro_Obstaculo - Posicion_Colision); #Vector que apunta de la posición "Posicion_Colision" al "Centro_Obstaculo"
	Normal_Colision = Normal_Colision/norm(Normal_Colision); #Normalizamos el vector previo
	Velocidad_Normal_Colision = dot(Normal_Colision, Velocidad_Colision)*Normal_Colision; #Componente normal a la recta tangente al obstáculo en el punto de colisión de la velocidad de la partícula al colisionar
	Velocidad_Colision = Velocidad_Colision - 2*Velocidad_Normal_Colision; #Velocidad de la partícula tras colisionar con el obstáculo
	Velocidad_Colision = Velocidad_Colision/norm(Velocidad_Colision); #Normalizamos la velocidad de la partícula tras la colisión 
	return Velocidad_Colision
end

#Función que calcula el tiempo de colisión, la posición de colisión y la velocidad tras la colisión de una partícula y un obstáculo.
#Posicion_Inicial: Arreglo con las coordenadas (Xi, Yi) asociadas a la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las coordenadas (Vxi, Vyi) asociadas a la velocidad inicial de la partícula.
#Radio_Particula: Radio de la circunferencia descrita por el desplazamiento de la partícula.
#Centro_Obstaculo: Arreglo con las coordenadas (Cx, Cy) asociadas al centro del obstáculo con el que la partícula podría colisionar.
#Radio_Obstaculo: Radio del obstáculo circular con el que la partícula podría colisionar
function colision_Condensada(Posicion_Inicial::Array, Velocidad_Inicial::Array, Radio_Particula::Float64, Centro_Obstaculo::Array, Radio_Obstaculo::Float64)
	Posicion_Colision, Tiempo_Colision = colision(Posicion_Inicial, Velocidad_Inicial, Radio_Particula, Centro_Obstaculo, Radio_Obstaculo); #Obtenemos la posición de la colisión y el tiempo para colisionar
	Centro_Particula = centro_Circunferencia(Velocidad_Inicial, Posicion_Inicial, Radio_Particula); #Obtenemos el centro de la circunferencia que describe la partícula al moverse
	Velocidad_Colision = velocidad_Particula(Posicion_Colision, Centro_Particula); #Obtenemos la velocidad de la partícula al colisionar
	Velocidad_Final = velocidad_Tras_Colision(Centro_Obstaculo, Posicion_Colision, Velocidad_Colision); #Obtenemos la velocidad de la partícula tras colisionar
	return Tiempo_Colision, Posicion_Colision, Velocidad_Final
end

#Función que nos regresa la posición de la partícula tras colisionar con un cilindro, su velocidad tras la colisión y el tiempo requerido para dicha colisión. Dependiendo de los parámetros que el usuario ingrese
#la función regresa también la información sobre el ángulo que forma la posición de la partícula al colisionar con el cilindro con respecto al centro del cilindro, el ángulo que forma la velocidad tras la colisión
#con respecto al centro del cilindro, la altura del plano en donde se dio la colisión y el centro del obstáculo en el plano 2D con el que colisionó la partícula.
#Posicion_Inicial: Arreglo con las coordenadas (Xi, Yi, Zi) asociadas a la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las coordenadas (Vxi, Vyi, Vzi) asociadas a la velocidad inicial de la partícula.
#Mat_Z_Cil: Matriz de rotación (0,0,1) -> EjeCil.
#Mat_Cil_Z: Matriz de rotación EjeCil -> (0,0,1).
#Centro_Cilindro: Arreglo con las coordenadas (Cx, Cy, Cz) asociadas a las coordenadas de la caja que contiene al cilindro de interés.
#Radio_Cilindro: Radio del cilindro de interés.
#Longitud_Cilindro: Longitud del cilindro de interés.
#Radio_Particula: Radio de la circunferencia descrita por el desplazamiento de la partícula.
#Poincare: Booleano que indica si se desea calcular la información necesaria para los mapeos de Poincaré (true) o no (false).
#Centros_Obstaculos: Booleano que indica si se desea calcular la información necesaria para graficar el obstáculo con el que colisiona la partícula (true) o no (false).
function colision_Cilindro(Posicion_Inicial::Array, Velocidad_Inicial::Array, Mat_Z_Cil::Array, Mat_Cil_Z::Array, Centro_Cilindro::Array, Radio_Cilindro::Float64, Longitud_Cilindro::Float64, Radio_Particula::Float64, Poincare::Bool, Centros_Obstaculos::Bool)
	Posicion_Cilindro_Vertical = Mat_Cil_Z*Posicion_Inicial; #Posición inicial de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
	Velocidad_Cilindro_Vertical = Mat_Cil_Z*Velocidad_Inicial; #Velocidad inicial de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
	Centro_Cilindro_Vertical = Mat_Cil_Z*Centro_Cilindro; #Coordenadas de la caja que contiene al cilindro de interés tras rotar el espacio "EjeCil" -> (0,0,1)
	Altura = Posicion_Cilindro_Vertical[3] - Centro_Cilindro_Vertical[3]; #Altura de la partícula en su posición inicial relativa al centro de la caja que contiene al cilindro tras rotar el espacio "EjeCil" -> (0,0,1)
	
	if abs(Altura) > Longitud_Cilindro/2 #Veamos si el plano de la trayectoria de la partícula corta al cilindro, si no es el caso se regresa la posicion inicial, la velocidad inicial y un tiempo infinito.
		if Poincare
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, Inf, Inf, Inf, Inf, [Inf, Inf] #Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision, θ_Pos_Col_Cen_Cil, θ_Velocidad, Altura, Centro_Cil_Vert_2D
			else
				return Posicion_Inicial, Velocidad_Inicial, Inf, Inf, Inf, Inf #Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision, θ_Pos_Col_Cen_Cil, θ_Velocidad, Altura
			end
		else
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, Inf, [Inf, Inf] #Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision, Centro_Cil_Vert_2D
			else
				return Posicion_Inicial, Velocidad_Inicial, Inf #Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision
			end
		end
	end

	Posicion_Cil_Vert_2D = [Posicion_Cilindro_Vertical[1], Posicion_Cilindro_Vertical[2]]; #Proyeccion al plano XY de la posicion de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
	Velocidad_Cil_Vert_2D = [Velocidad_Cilindro_Vertical[1], Velocidad_Cilindro_Vertical[2]]; #Proyeccion al plano XY de la velocidad de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
	Velocidad_Cil_Vert_2D = Velocidad_Cil_Vert_2D/norm(Velocidad_Cil_Vert_2D); #Normalizamos el vector anterior
	Centro_Cil_Vert_2D = [Centro_Cilindro_Vertical[1], Centro_Cilindro_Vertical[2]]; #Proyección al plano XY de las coordenadas de la caja que contiene al cilindro tras rotar el espacio "EjeCil" -> (0,0,1)
	Tiempo_Colision, Posicion_Colision, Velocidad_Colision = colision_Condensada(Posicion_Cil_Vert_2D, Velocidad_Cil_Vert_2D, Radio_Particula, Centro_Cil_Vert_2D, Radio_Cilindro);
	Posicion_Colision_Centrada_Cilindro = Posicion_Colision - Centro_Cil_Vert_2D; #Obtenemos la posición de la colisión de la partícula relativa al centro del obstáculo con el que colisionó

	if Poincare
		θ_Pos_Col_Cen_Cil = atan(Posicion_Colision_Centrada_Cilindro[2], Posicion_Colision_Centrada_Cilindro[1]); #Angulo con respecto al eje X de la posición de colisión relativa al centro del obstáculo
		θ_Velocidad = atan(Velocidad_Colision[2], Velocidad_Colision[1]); #Angulo con respecto al eje X de la velocidad tras la colisión relativa al centro del obstáculo
	end

	Posicion_Colision_3D = push!(Posicion_Colision, Posicion_Cilindro_Vertical[3]); #Agregamos la altura del plano en el que está contenida la partícula a la posición 2D en dicho plano de la partícula
	Velocidad_Colision_3D = push!(Velocidad_Colision, Velocidad_Cilindro_Vertical[3]); #Agregamos la altura del plano en el que está contenida la partícula a la velocidad 2D en dicho plano de la partícula
	Posicion_Colision_3D = Mat_Z_Cil*Posicion_Colision_3D; #Obtenemos la posición de la partícula al colisionar con el cilindro tras rotar el espacio (0,0,1) -> "EjeCil"
	Velocidad_Colision_3D = Mat_Z_Cil*Velocidad_Colision_3D; #Obtenemos la velocidad de la partícula tras colisionar con el cilindro tras rotar el espacio (0,0,1) -> "EjeCil"
	Velocidad_Colision_3D = Velocidad_Colision_3D/norm(Velocidad_Colision_3D); #Normalizamos la velocidad previa

	if Poincare
		if Centros_Obstaculos
			return Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision, θ_Pos_Col_Cen_Cil, θ_Velocidad, Altura, Centro_Cil_Vert_2D 
		else
			return Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision, θ_Pos_Col_Cen_Cil, θ_Velocidad, Altura
		end
	else
		if Centros_Obstaculos
			return Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision, Centro_Cil_Vert_2D
		else
			return Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision
		end
	end
end

#Función que nos regresa la posición y la velocidad de la partícula tras un tiempo de vuelo Δt.
#Posicion_Inicial: Arreglo con las coordenadas (Xi, Yi, Zi) asociadas a la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las coordenadas (Vxi, Vyi, Vzi) asociadas a la velocidad inicial de la partícula.
#Mat_Z_Cil: Matriz de rotación (0,0,1) -> EjeCil.
#Mat_Cil_Z: Matriz de rotación EjeCil -> (0,0,1).
#Radio_Particula: Radio de la circunferencia descrita por el desplazamiento de la partícula.
#Δt: Tiempo de vuelo de la partícula.
function avanza(Posicion_Inicial::Array, Velocidad_Inicial::Array, Mat_Z_Cil::Array, Mat_Cil_Z::Array, Radio_Particula::Float64, Δt::Float64)
	Posicion_Cilindro_Vertical = Mat_Cil_Z*Posicion_Inicial; #Posición inicial de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
	Velocidad_Cilindro_Vertical = Mat_Cil_Z*Velocidad_Inicial; #Velocidad inicial de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)

	Posicion_Cil_Vert_2D = [Posicion_Cilindro_Vertical[1], Posicion_Cilindro_Vertical[2]]; #Proyección al plano XY de la posición inicial de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
	Velocidad_Cil_Vert_2D = [Velocidad_Cilindro_Vertical[1], Velocidad_Cilindro_Vertical[2]]; #Proyección al plano XY de la velocidad inicial de la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
	Velocidad_Cil_Vert_2D = Velocidad_Cil_Vert_2D/norm(Velocidad_Cil_Vert_2D); #Normalizamos el vector previo

	Posicion_Intermedia, Velocidad_Intermedia = avance(Posicion_Cil_Vert_2D, Velocidad_Cil_Vert_2D, Radio_Particula, Δt); #Obtenemos la posición y la velocidad de la partícula tras desplazarse un tiempo Δt

	push!(Posicion_Intermedia, Posicion_Cilindro_Vertical[3]); #Agregamos al vector de la posición intermedia de la partícula la altura del plano donde se mueve la partícula
	push!(Velocidad_Intermedia, Velocidad_Cilindro_Vertical[3]); #Agremamos al vector de la velocidad intermedia de la partícula la componente en Z de su velocidad

	Posicion_Final = Mat_Z_Cil*Posicion_Intermedia; #Posición final de la partícula tras avanzar un tiempo Δt y rotar el espacio (0,0,1) -> "EjeCil"
	Velocidad_Final = Mat_Z_Cil*Velocidad_Intermedia; #Velocidad final de la partícula tras avanzar un tiempo Δt y rotar el espacio (0,0,1) -> "EjeCil"
	Velocidad_Final = Velocidad_Final/norm(Velocidad_Final); #Normalizamos el vector previo

	return Posicion_Final, Velocidad_Final
end

#Función que regresa la posición final de la partícula, su velocidad final y las coordenadas (X,Y,Z) de la caja que contendría a la partícula si no hubiera condiciones periódicas a la frontera. Dependiendo de los
#parámetros del usuario, puede regresar de igual forma los ángulos con respecto al eje X de la colisión de la partícula con respecto al centro de cada uno de los obstáculos con los que colisionó así como los ángulos 
#con respecto al eje X de la velocidad tras la colisión con cada uno de los obstaculos con los que colisionó (Poincare = true); de igual forma puede regresar los centros de los obstáculos en el plano XY con los que
#colisionó (Centros_Obstaculos = true).
#Posicion_Inicial: Arreglo con las coordenadas (Xi, Yi, Zi) asociadas a la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las coordenadas (Vxi, Vyi, Vzi) asociadas a la velocidad inicial de la partícula.
#Mat_Z_Cil: Matriz de rotación (0,0,1) -> EjeCil.
#Mat_Cil_Z: Matriz de rotación EjeCil -> (0,0,1).
#LongitudCil: Longitud del cilindro que colocamos dentro de la caja central.
#LadoCaja: Longitud de la arista de la caja.
#Radio_Particula: Radio de la circunferencia descrita por el desplazamiento de la partícula.
#Radio_Cilindro: Radio del cilindro que colocamos dentro de la caja central.
#Tiempo_Vuelo: Tiempo de vuelo total que tendrá la partícula.
#Centros: Arreglo en donde guardaremos los centros de los obstáculos con los que la partícula va a colisionar (SÓLO TIENE UTILIDAD SI Centros_Obstaculos = true).
#N_Caja: Coordenadas (X,Y,Z) de la caja que contendría a la partícula en su posición inicial si no hubiera condiciones periódicas a la frontera.
#Poincare: Booleano que determina si queremos calcular o no los parámetros relacionados al mapeo de Poincaré.
#Centros_Obstaculos: Booleano que determina si queremos calcular o no los parámetros relacionados a los centros en el plano de los obstáculos con los que la partícula va a colisionar.
#Eficiencia: Booleano para determinar si hacemos los cálculos considerando 3 cajas o las 27 cajas.
function lorentz_Cuasi_Magnetico(Posicion_Inicial::Array, Velocidad_Inicial::Array, Mat_Z_Cil::Array, Mat_Cil_Z::Array, LongitudCil::Float64, LadoCaja::Float64, Radio_Particula::Float64, Radio_Cilindro::Float64, Tiempo_Vuelo::Float64, Centros::Array, N_Caja::Array, Poincare::Bool, Centros_Obstaculos::Bool, Eficiencia::Bool)
	if Eficiencia #Sólo se considera la caja central y la caja arriba y debajo de la caja central en el cálculo de posibles colisiones
		N_Total = zeros(3); #Arreglo que contendrá las coordenadas (Cx,Cy,Cz) de la caja que contendría a la partícula si no hubiera condiciones periódicas a la frontera
		δT = 0.1; #Tiempo vuelo discreto (Para graficar trayectorias)
		i = 0;

		if Poincare
			Θc = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la posición de colisión de la partícula con respecto al centro del obstáculo
			Θv = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la velocidad tras la colisión de la partícula con respecto al centro del obstáculo
		end

		while Tiempo_Vuelo > 0.0
			i += 1;
			t = Inf; #Parámetro que funcionará como "el tiempo que tarda en realizarse la próxima colisión de la partícula"
			Posicion_Inicial, Coordenada_Caja = caja_Periodica!(Posicion_Inicial, LadoCaja); #Obtenemos las posición relativa a la caja central de la partícula, así como las coordenadas de la caja que la contendría
			N_Total += Coordenada_Caja; #Sumamos a las coordenadas de la caja que contendrá a la partícula, las coordenadas (relativas a la caja central) de la caja que la contiene en este paso

			Posicion_1 = copy(Posicion_Inicial); #Hacemos una copia de la posición inicial de la partícula
			Velocidad_1 = copy(Velocidad_Inicial); #Hacemos una copia de la velocidad inicial de la partícula
			
			if Centros_Obstaculos
				Centros_Cilindros = [0.0, 0.0]; #Arreglo que contendrá los centros de los obstáculos con los que la partícula colisiona
			end

			Coordenada_X_Cil = 0; #Coordenada X de la caja donde se encuentra el cilindro con el que queremos verificar si colisiona la partícula o no
			Coordenada_Y_Cil = 0; #Coordenada Y de la caja donde se encuentra el cilindro con el que queremos verificar si colisiona la partícula o no

			for Coordenada_Z_Cil in -1:1
				if Poincare
					if Centros_Obstaculos
						Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					else
						Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					end
				else
					if Centros_Obstaculos
						Posicion_2, Velocidad_2, Tiempo_Colision_1, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					else
						Posicion_2, Velocidad_2, Tiempo_Colision_1 = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					end
				end

				if Tiempo_Colision_1 < t #Calcula la colisión que consuma menos tiempo de vuelo (la primera colisión)
					t = Tiempo_Colision_1;
					Posicion_1 = copy(Posicion_2);
					Velocidad_1 = copy(Velocidad_2);

					if Centros_Obstaculos
						Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);
					end

					if Poincare
						push!(Θc, θ_Pos_Col_Cen_Cil_1);
						push!(Θv, θ_Velocidad_1);
					end
				end
			end

			if t < δT #Si la próxima colisión se da en un tiempo menor que el parámetro de tiempo discreto
				Tiempo_Vuelo -= t; #Réstale al tiempo de vuelo el tiempo que tarda en darse la próxima colisión
				Posicion_Inicial = Posicion_1; #Actualiza la posición inicial de la partícula a la posición de la próxima colisión
				Velocidad_Inicial = Velocidad_1; #Actualiza la velocidad inicial de la partícula a la velocidad tras la próxima colisión
				
				if Centros_Obstaculos
					NN = Mat_Cil_Z*(N_Total + N_Caja); #Obtiene las coordenadas de la caja final que contiene a la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
					push!(Centros, Centros_Cilindros + [NN[1], NN[2]]); #Agrega al arreglo Centros las coordenadas (X,Y) del obstáculo con el que la partícula colisionó si no hubiera condiciones periódicas
				end
			else #Si la colisión se da después del parámetro de tiempo discreto
				Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, Radio_Particula, δT); #Evolucionamos la trayectoria de la partícula la unidad de tiempo discreto
				Tiempo_Vuelo -= δT; #Le quitamos al tiempo de vuelo el tiempo invertido en desplazar la partícula de manera discreta
			end
		end

		if Poincare
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv
			end
		else
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total
			end
		end
	else
		N_Total = zeros(3); #Arreglo que contendrá las coordenadas (Cx,Cy,Cz) de la caja que contendría a la partícula si no hubiera condiciones periódicas a la frontera
		δT = 0.1; #Tiempo vuelo discreto (Para graficar trayectorias)
		i = 0;

		if Poincare
			Θc = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la posición de colisión de la partícula con respecto al centro del obstáculo
			Θv = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la velocidad tras la colisión de la partícula con respecto al centro del obstáculo
		end

		while Tiempo_Vuelo > 0.0
			i += 1;
			t = Inf; #Parámetro que funcionará como "el tiempo que tarda en realizarse la próxima colisión de la partícula"
			Posicion_Inicial, Coordenada_Caja = caja_Periodica!(Posicion_Inicial, LadoCaja); #Obtenemos las posición relativa a la caja central de la partícula, así como las coordenadas de la caja que la contendría
			N_Total += Coordenada_Caja; #Sumamos a las coordenadas de la caja que contendrá a la partícula, las coordenadas (relativas a la caja central) de la caja que la contiene en este paso

			Posicion_1 = copy(Posicion_Inicial); #Hacemos una copia de la posición inicial de la partícula
			Velocidad_1 = copy(Velocidad_Inicial); #Hacemos una copia de la velocidad inicial de la partícula

			if Centros_Obstaculos
				Centros_Cilindros = [0.0, 0.0]; #Arreglo que contendrá los centros de los obstáculos con los que la partícula colisiona
			end

			for ix in -1:1
				for iy in -1:1
					for iz in -1:1
						if Poincare
							if Centros_Obstaculos
								Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							else
								Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							end
						else
							if Centros_Obstaculos
								Posicion_2, Velocidad_2, Tiempo_Colision_1, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							else
								Posicion_2, Velocidad_2, Tiempo_Colision_1 = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							end
						end
                        
						if Tiempo_Colision_1 < t

							t = Tiempo_Colision_1;
							Posicion_1 = copy(Posicion_2);
							Velocidad_1 = copy(Velocidad_2);

							if Centros_Obstaculos
								Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);
							end

							if Poincare
								push!(Θc, θ_Pos_Col_Cen_Cil_1);
								push!(Θv, θ_Velocidad_1);
							end
						end
					end
				end
			end

			if t < δT #Si la próxima colisión se da en un tiempo menor que el parámetro de tiempo discreto
				Tiempo_Vuelo -= t; #Réstale al tiempo de vuelo el tiempo que tarda en darse la próxima colisión
				Posicion_Inicial = Posicion_1; #Actualiza la posición inicial de la partícula a la posición de la próxima colisión
				Velocidad_Inicial = Velocidad_1; #Actualiza la velocidad inicial de la partícula a la velocidad tras la próxima colisión
				
				if Centros_Obstaculos
					NN = Mat_Cil_Z*(N_Total + N_Caja); #Obtiene las coordenadas de la caja final que contiene a la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
					push!(Centros, Centros_Cilindros + [NN[1], NN[2]]); #Agrega al arreglo Centros las coordenadas (X,Y) del obstáculo con el que la partícula colisionó si no hubiera condiciones periódicas
				end
			else #Si la colisión se da después del parámetro de tiempo discreto
				Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, Radio_Particula, δT); #Evolucionamos la trayectoria de la partícula la unidad de tiempo discreto
				Tiempo_Vuelo -= δT; #Le quitamos al tiempo de vuelo el tiempo invertido en desplazar la partícula de manera discreta
			end
		end

		if Poincare
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv
			end
		else
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total
			end
		end
	end
end

#Función que regresa la posición final de la partícula, su velocidad final y las coordenadas (X,Y,Z) de la caja que contendría a la partícula si no hubiera condiciones periódicas a la frontera. Dependiendo de los
#parámetros del usuario, puede regresar de igual forma los ángulos con respecto al eje X de la colisión de la partícula con respecto al centro de cada uno de los obstáculos con los que colisionó así como los ángulos 
#con respecto al eje X de la velocidad tras la colisión con cada uno de los obstaculos con los que colisionó (Poincare = true); de igual forma puede regresar los centros de los obstáculos en el plano XY con los que
#colisionó (Centros_Obstaculos = true).
#Posicion_Inicial: Arreglo con las coordenadas (Xi, Yi, Zi) asociadas a la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las coordenadas (Vxi, Vyi, Vzi) asociadas a la velocidad inicial de la partícula.
#Mat_Z_Cil: Matriz de rotación (0,0,1) -> EjeCil.
#Mat_Cil_Z: Matriz de rotación EjeCil -> (0,0,1).
#LongitudCil: Longitud del cilindro que colocamos dentro de la caja central.
#LadoCaja: Longitud de la arista de la caja.
#Radio_Particula: Radio de la circunferencia descrita por el desplazamiento de la partícula.
#Radio_Cilindro: Radio del cilindro que colocamos dentro de la caja central.
#Tiempo_Vuelo: Tiempo de vuelo total que tendrá la partícula.
#Centros: Arreglo en donde guardaremos los centros de los obstáculos con los que la partícula va a colisionar (SÓLO TIENE UTILIDAD SI Centros_Obstaculos = true).
#N_Caja: Coordenadas (X,Y,Z) de la caja que contendría a la partícula en su posición inicial si no hubiera condiciones periódicas a la frontera.
#Poincare: Booleano que determina si queremos calcular o no los parámetros relacionados al mapeo de Poincaré.
#Centros_Obstaculos: Booleano que determina si queremos calcular o no los parámetros relacionados a los centros en el plano de los obstáculos con los que la partícula va a colisionar.
#Eficiencia: Booleano para determinar si hacemos los cálculos considerando 3 cajas o las 27 cajas.
function lorentz_Cuasi_Magnetico2(Posicion_Inicial::Array, Velocidad_Inicial::Array, Mat_Z_Cil::Array, Mat_Cil_Z::Array, LongitudCil::Float64, LadoCaja::Float64, Radio_Particula::Float64, Radio_Cilindro::Float64, Tiempo_Vuelo::Float64, Centros::Array, N_Caja::Array, Poincare::Bool, Centros_Obstaculos::Bool, Eficiencia::Bool)
    δ = min(.5,Tiempo_Vuelo)
    δ2 = min(.1, Tiempo_Vuelo)
	if Eficiencia #Sólo se considera la caja central y la caja arriba y debajo de la caja central en el cálculo de posibles colisiones
		N_Total = zeros(3); #Arreglo que contendrá las coordenadas (Cx,Cy,Cz) de la caja que contendría a la partícula si no hubiera condiciones periódicas a la frontera
		δT = min(1000.1, Tiempo_Vuelo);
		i = 0;

		if Poincare
			Θc = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la posición de colisión de la partícula con respecto al centro del obstáculo
			Θv = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la velocidad tras la colisión de la partícula con respecto al centro del obstáculo
			H = []; #Arreglo donde irá la altura donde la partícula colisiona con el cilindro respecto al centro del cilindro.
		end

		while Tiempo_Vuelo > 0.0
			i += 1;
			t = Inf; #Parámetro que funcionará como "el tiempo que tarda en realizarse la próxima colisión de la partícula"
			Posicion_Inicial, Coordenada_Caja = caja_Periodica!(Posicion_Inicial, LadoCaja); #Obtenemos las posición relativa a la caja central de la partícula, así como las coordenadas de la caja que la contendría
			N_Total += Coordenada_Caja; #Sumamos a las coordenadas de la caja que contendrá a la partícula, las coordenadas (relativas a la caja central) de la caja que la contiene en este paso

			Posicion_1 = copy(Posicion_Inicial); #Hacemos una copia de la posición inicial de la partícula
			Velocidad_1 = copy(Velocidad_Inicial); #Hacemos una copia de la velocidad inicial de la partícula
			
			if Centros_Obstaculos
				Centros_Cilindros = [0.0, 0.0]; #Arreglo que contendrá los centros de los obstáculos con los que la partícula colisiona
			end

			Coordenada_X_Cil = 0; #Coordenada X de la caja donde se encuentra el cilindro con el que queremos verificar si colisiona la partícula o no
			Coordenada_Y_Cil = 0; #Coordenada Y de la caja donde se encuentra el cilindro con el que queremos verificar si colisiona la partícula o no

			if Poincare
            	θ_Pos_Col_Cen_Cil_2 = 0;
            	θ_Velocidad_2 = 0;
				h2 = 0
			end

			for Coordenada_Z_Cil in [-1,0,1]
				if Poincare
					if Centros_Obstaculos
						Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					else
						Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					end
				else
					if Centros_Obstaculos
						Posicion_2, Velocidad_2, Tiempo_Colision_1, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					else
						Posicion_2, Velocidad_2, Tiempo_Colision_1 = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
					end
				end
				
				if Tiempo_Colision_1 < t #Calcula la colisión que consuma menos tiempo de vuelo (la primera colisión)
					t = Tiempo_Colision_1;
					Posicion_1 = copy(Posicion_2);
					Velocidad_1 = copy(Velocidad_2);

					if Centros_Obstaculos
						Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);
					end

					if Poincare
						θ_Velocidad_2 = copy(θ_Velocidad_1);
						θ_Pos_Col_Cen_Cil_2 = copy(θ_Pos_Col_Cen_Cil_1);
						h2 = copy(Altura);
					end
				end
			end
			
			if Poincare
				push!(H,h2)
				push!(Θc, θ_Pos_Col_Cen_Cil_2);
				push!(Θv, θ_Velocidad_2);
			end

			if t < δT #Si la próxima colisión se da en un tiempo menor que el parámetro de tiempo discreto
				Tiempo_Vuelo -= t; #Réstale al tiempo de vuelo el tiempo que tarda en darse la próxima colisión
				Posicion_Inicial = Posicion_1; #Actualiza la posición inicial de la partícula a la posición de la próxima colisión
				Velocidad_Inicial = Velocidad_1; #Actualiza la velocidad inicial de la partícula a la velocidad tras la próxima colisión
				
				if Centros_Obstaculos
					NN = Mat_Cil_Z*(N_Total + N_Caja); #Obtiene las coordenadas de la caja final que contiene a la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
					push!(Centros, Centros_Cilindros + [NN[1], NN[2]]); #Agrega al arreglo Centros las coordenadas (X,Y) del obstáculo con el que la partícula colisionó si no hubiera condiciones periódicas
				end
			else #Si la colisión se da después del parámetro de tiempo discreto
				Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, Radio_Particula, δ2) #Evolucionamos la trayectoria de la partícula la unidad de tiempo discreto
				Tiempo_Vuelo -= δ2; #Le quitamos al tiempo de vuelo el tiempo invertido en desplazar la partícula de manera discreta
			end
		end

		if Poincare
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, H, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, H
			end
		else
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total
			end
		end
	else
		N_Total = zeros(3); #Arreglo que contendrá las coordenadas (Cx,Cy,Cz) de la caja que contendría a la partícula si no hubiera condiciones periódicas a la frontera
		δT = min(1000.1, Tiempo_Vuelo);
		i = 0;

		if Poincare
			Θc = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la posición de colisión de la partícula con respecto al centro del obstáculo
			Θv = []; #Arreglo donde irá el ángulo que forma con respecto al eje X la velocidad tras la colisión de la partícula con respecto al centro del obstáculo
			H = []; #Arreglo donde irá la altura donde la partícula colisiona con el cilindro respecto al centro del cilindro.
		end

		while Tiempo_Vuelo > 0.0
			i += 1;
			t = Inf; #Parámetro que funcionará como "el tiempo que tarda en realizarse la próxima colisión de la partícula"
			Posicion_Inicial, Coordenada_Caja = caja_Periodica!(Posicion_Inicial, LadoCaja); #Obtenemos las posición relativa a la caja central de la partícula, así como las coordenadas de la caja que la contendría
			N_Total += Coordenada_Caja; #Sumamos a las coordenadas de la caja que contendrá a la partícula, las coordenadas (relativas a la caja central) de la caja que la contiene en este paso

			Posicion_1 = copy(Posicion_Inicial); #Hacemos una copia de la posición inicial de la partícula
			Velocidad_1 = copy(Velocidad_Inicial); #Hacemos una copia de la velocidad inicial de la partícula
			
			if Centros_Obstaculos
				Centros_Cilindros = [0.0, 0.0]; #Arreglo que contendrá los centros de los obstáculos con los que la partícula colisiona
			end

			if Poincare
            	θ_Pos_Col_Cen_Cil_2 = 0;
            	θ_Velocidad_2 = 0;
				h2 = 0
			end

			for ix in -1:1
				for iy in -1:1
					for iz in -1:1
						if Poincare
							if Centros_Obstaculos
								Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							else
								Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, Altura = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							end
						else
							if Centros_Obstaculos
								Posicion_2, Velocidad_2, Tiempo_Colision_1, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							else
								Posicion_2, Velocidad_2, Tiempo_Colision_1 = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, [ix,iy,iz], Radio_Cilindro, LongitudCil, Radio_Particula, Poincare, Centros_Obstaculos);
							end
						end
						
						if Tiempo_Colision_1 < t #Calcula la colisión que consuma menos tiempo de vuelo (la primera colisión)
							t = Tiempo_Colision_1;
							Posicion_1 = copy(Posicion_2);
							Velocidad_1 = copy(Velocidad_2);

							if Centros_Obstaculos
								Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);
							end

                            if Poincare
								θ_Velocidad_2 = copy(θ_Velocidad_1);
								θ_Pos_Col_Cen_Cil_2 = copy(θ_Pos_Col_Cen_Cil_1);
								h2 = copy(Altura);
							end
						end
					end
				end
			end

			if Poincare
				push!(H,h2)
				push!(Θc, θ_Pos_Col_Cen_Cil_2);
				push!(Θv, θ_Velocidad_2);
			end

			if t < δT #Si la próxima colisión se da en un tiempo menor que el parámetro de tiempo discreto
				Tiempo_Vuelo -= t; #Réstale al tiempo de vuelo el tiempo que tarda en darse la próxima colisión
				Posicion_Inicial = Posicion_1; #Actualiza la posición inicial de la partícula a la posición de la próxima colisión
				Velocidad_Inicial = Velocidad_1; #Actualiza la velocidad inicial de la partícula a la velocidad tras la próxima colisión
				
				if Centros_Obstaculos
					NN = Mat_Cil_Z*(N_Total + N_Caja); #Obtiene las coordenadas de la caja final que contiene a la partícula tras rotar el espacio "EjeCil" -> (0,0,1)
					push!(Centros, Centros_Cilindros + [NN[1], NN[2]]); #Agrega al arreglo Centros las coordenadas (X,Y) del obstáculo con el que la partícula colisionó si no hubiera condiciones periódicas
				end
			else #Si la colisión se da después del parámetro de tiempo discreto
                Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, Mat_Z_Cil, Mat_Cil_Z, Radio_Particula, δ)
				Tiempo_Vuelo -= δ;
			end
		end

		if Poincare
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, H, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, H
			end
		else
			if Centros_Obstaculos
				return Posicion_Inicial, Velocidad_Inicial, N_Total, Centros
			else
				return Posicion_Inicial, Velocidad_Inicial, N_Total
			end
		end
	end
end

#Función que nos regresa la gráfica de un círculo dado su radio y su centro
function dibuja_Circulo(Radio_Circulo::Float64, Centro_Circulo::Array, color = "red")
	X = [Radio_Circulo];
	Y = [0.0];

	for i in 1:99
		x = Radio_Circulo*cos(2*π/100 * i);
		y = Radio_Circulo*sin(2*π/100 * i);
		push!(X, x);
		push!(Y, y);
	end

	push!(X, Radio_Circulo);
	push!(Y, 0.0);

	X .+= Centro_Circulo[1];
	Y .+= Centro_Circulo[2];

	return plot!(X, Y, color = color)
end
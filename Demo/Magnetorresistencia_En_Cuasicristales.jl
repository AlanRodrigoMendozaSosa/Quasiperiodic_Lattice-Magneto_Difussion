using LinearAlgebra #A partir de la versión 1.0 de Julia, las funciones "norm" y "dot" se encuentran dentro de dicha paquetería.

function big_I_3x3()
	reshape([1,0,0,0,1,0,0,0,1],3,3)
end

function skew_Symm(A)
	reshape([0, A[3], -A[2], -A[3], 0, A[1], A[2], -A[1], 0],3,3);
end

function caja_Periodica(Posicion, LadoCaja) #El vector posición se verá modificado al finalizar la función.
	#La caja está centrada en el origen del sistema de referencia y tiene por semilado LadoCaja/2.
	N = length(Posicion);
	Caja = floor.((Posicion .+ LadoCaja/2.) ./ LadoCaja) * LadoCaja;#Nos indica la coordenada de la caja en la cual se encuentra la partícula antes de las condiciones periódicas.

	#Condición de que la partícula se mantenga dentro de la caja en cada coordenada.
	for i in 1:N
		if (-LadoCaja/2. > Posicion[i]) || (Posicion[i] >= LadoCaja/2.)
			Posicion[i] = mod(Posicion[i] + LadoCaja/2., LadoCaja) - LadoCaja/2.;
		end
	end

	return Posicion, Caja
end

function mat_Rot(Normal) #Matriz que transforma al vector Z en el vector Normal.
	Tipo = typeof(Normal[1]);
	E1 = [parse(Tipo, "0"), parse(Tipo, "0"), parse(Tipo, "1")]; #vector unitario en Z.
	E2 = Normal/norm(Normal); #vector unitario en dirección del eje del cilindro.

	Cruz = cross(E1,E2);
	Coseno = dot(E1,E2);
	(Coseno == -1) ? throw(ArgumentError("La normal no puede ser colineal al eje Z")) : nothing;
	Discriminante = 1/(1 + Coseno);

	Matriz_Rotacion = big_I_3x3() + skew_Symm(Cruz) + Discriminante*(skew_Symm(Cruz)*skew_Symm(Cruz));
end

function vel_Rand(Matriz_Rotacion)
	θ = (2*π)*rand();
	Vel_Sin_Girar = [cos(θ), sin(θ), 0];
	Matriz_Rotacion*Vel_Sin_Girar;
end

function centros_Cajas_3D(LadoCaja)
	Tipo = typeof(LadoCaja);
	A = [[-LadoCaja, LadoCaja, -LadoCaja], 
        [Tipo(0), LadoCaja, -LadoCaja], 
        [LadoCaja, LadoCaja, -LadoCaja], 
        [-LadoCaja, Tipo(0), -LadoCaja], 
        [Tipo(0), Tipo(0), -LadoCaja], 
        [LadoCaja, Tipo(0), -LadoCaja], 
        [-LadoCaja, -LadoCaja, -LadoCaja], 
        [Tipo(0), -LadoCaja, -LadoCaja], 
        [LadoCaja, -LadoCaja, -LadoCaja], 
        [-LadoCaja, LadoCaja, Tipo(0)], 
        [Tipo(0), LadoCaja, Tipo(0)], 
        [LadoCaja, LadoCaja, Tipo(0)], 
        [-LadoCaja, Tipo(0), Tipo(0)], 
        [Tipo(0),Tipo(0),Tipo(0)], 
        [LadoCaja, Tipo(0), Tipo(0)], 
        [-LadoCaja, -LadoCaja, Tipo(0)], 
        [Tipo(0), -LadoCaja, Tipo(0)], 
        [LadoCaja, -LadoCaja, Tipo(0)], 
        [-LadoCaja, LadoCaja, LadoCaja], 
        [Tipo(0), LadoCaja, LadoCaja], 
        [LadoCaja, LadoCaja, LadoCaja], 
        [-LadoCaja, Tipo(0), LadoCaja], 
        [Tipo(0), Tipo(0), LadoCaja], 
        [LadoCaja, Tipo(0), LadoCaja], 
        [-LadoCaja, -LadoCaja, LadoCaja], 
        [Tipo(0), -LadoCaja, LadoCaja], 
        [LadoCaja, -LadoCaja, LadoCaja]]
    return A
end

function longitud_Cilindro(LadoCaja, Matriz_Rotacion_Inv)
	#Definimos las coordenadas de las esquinas del cubo sin rotar.
    Arreglo_Aristas_No_Rot = [[-(LadoCaja)/2, (LadoCaja)/2, -(LadoCaja)/2], 
        [-(LadoCaja)/2, -(LadoCaja)/2, -(LadoCaja)/2], 
        [(LadoCaja)/2, (LadoCaja)/2, -(LadoCaja)/2], 
        [(LadoCaja)/2, -(LadoCaja)/2, -(LadoCaja)/2],
        [-(LadoCaja)/2, (LadoCaja)/2, (LadoCaja)/2], 
        [-(LadoCaja)/2, -(LadoCaja)/2, (LadoCaja)/2], 
        [(LadoCaja)/2, (LadoCaja)/2, (LadoCaja)/2], 
        [(LadoCaja)/2, -(LadoCaja)/2, (LadoCaja)/2]];

    Tipo = typeof(LadoCaja)

    #Giramos el cubo de modo que el cilindro quede alineado con el eje Z, esto nos da las siguientes coordenadas.
    Arreglo_Aristas_Rot = [];

    for i in 1:length(Arreglo_Aristas_No_Rot)
    	x = Matriz_Rotacion_Inv*Arreglo_Aristas_No_Rot[i];
    	push!(Arreglo_Aristas_Rot, x);
    end

    #Obtenemos la máxima longitud entre cada uno de los puntos en su coordenada en Z
    Arreglo_Longitud_Z = Tipo[];

    for i in 1:(length(Arreglo_Aristas_Rot)-1)
    	for j in (i+1):length(Arreglo_Aristas_Rot)
    		LongitudZ = abs(Arreglo_Aristas_Rot[i][3] - Arreglo_Aristas_Rot[j][3]);
    		push!(Arreglo_Longitud_Z, LongitudZ);
    	end
    end

    Maxima_Longitud_Z = maximum(Arreglo_Longitud_Z);

    return Maxima_Longitud_Z
end

function posiciones_Iniciales(EjeCil, RadioCil, LadoCaja, Iteraciones)
	#EjeCil es el vector que apunta al eje del cilindro: DEBE ESTAR EN EL SEMIEJE SUPERIOR (EjeCil.z > 0)
	#LadoCaja es la longitud de un lado completo de la caja.

	Eje_Cil_Uni = EjeCil/norm(EjeCil); #Vector unitario y paralelo al eje del cilindro.
	Mat_Z_Cil = mat_Rot(EjeCil); #Matriz de rotación del eje Z al eje del cilindro.
	Mat_Cil_Z = (Mat_Z_Cil)^-1; #Matriz de rotación del eje del cilindro al eje Z.
	LadoCil = longitud_Cilindro(LadoCaja, Mat_Cil_Z); #Longitud del cilindro.
	Tipo = typeof(LadoCaja);

	#Definiendo los centros de todas las cajas adyacentes a la caja central.
	Arreglo_Centros_Cilindros = centros_Cajas_3D(LadoCaja);

	#Obtengamos estos centros tras rotar el espacio y tener los cilindros verticales.
	Arreglo_Centros_Cilindros_Rotados = [];

	for j in 1:length(Arreglo_Centros_Cilindros)
		x = Mat_Cil_Z*Arreglo_Centros_Cilindros[j];
		push!(Arreglo_Centros_Cilindros_Rotados, x);
	end

	Arreglo_Posiciones = []; #Arreglo donde irán las coordenadas de puntos válidos.

	for i in 1:Iteraciones
		#Llaves para iniciar el código.
		Calcular = true;
		Intentos = 0;

		try
			while Calcular == true && Intentos < 1e4
				Intentos += 1;

				#En este punto coloquemos una partícula, sin importarnos si estamos dentro o fuera del cilindro.
				PrevX = Tipo(rand() - parse(Tipo, "0.5"));
				PrevY = Tipo(rand() - parse(Tipo, "0.5"));
				PrevZ = Tipo(rand() - parse(Tipo, "0.5"));
				Prev_Pos = [PrevX, PrevY, PrevZ];

				#Calculemos si ese punto es un punto válido.
				#Primero rotemos para que los cilindros queden verticales.
				Pos = Mat_Cil_Z*Prev_Pos;

				Llaves = Bool[];

				#Ahora coloquemos el origen en el centro de cada cilindro rotado y veamos si el punto cae dentro.
				for j in 1:length(Arreglo_Centros_Cilindros_Rotados)
					Posicion_Punto = Pos - Arreglo_Centros_Cilindros_Rotados[j]; #Punto de la partícula a prueba tras la rotación y traslación espacial del cilindro J al origen.
					((Posicion_Punto[3] < LadoCil/2 && Posicion_Punto[3] > -(LadoCil/2)) && ((Posicion_Punto[1])^2 + (Posicion_Punto[2])^2) < (RadioCil)^2) ? (push!(Llaves, true)) : (push!(Llaves, false));
				end

				any(Llaves) ? (nothing) : (push!(Arreglo_Posiciones, Prev_Pos); Calcular = false);

				Intentos == 1e4 ? throw("Demasiados Intentos sin encontrar un punto válido") : nothing;
			end
		catch
			nothing
		end
	end

	return Arreglo_Posiciones
end

function posiciones_Iniciales_2(EjeCil, RadioCil, LadoCaja, Iteraciones)
    #EjeCil es el vector que apunta al eje del cilindro: DEBE ESTAR EN EL SEMIEJE SUPERIOR (EjeCil.z > 0)
	#LadoCaja es la longitud de un lado completo de la caja.

	Eje_Cil_Uni = EjeCil/norm(EjeCil); #Vector unitario y paralelo al eje del cilindro.
	Mat_Z_Cil = mat_Rot(EjeCil); #Matriz de rotación del eje Z al eje del cilindro.
	Mat_Cil_Z = (Mat_Z_Cil)^-1; #Matriz de rotación del eje del cilindro al eje Z.
	LadoCil = longitud_Cilindro(LadoCaja, Mat_Cil_Z); #Longitud del cilindro.
	Tipo = typeof(LadoCaja);

	#Definiendo los centros de todas las cajas adyacentes a la caja central.
	Arreglo_Centros_Cilindros = centros_Cajas_3D(LadoCaja);

	#Obtengamos estos centros tras rotar el espacio y tener los cilindros verticales.
	Arreglo_Centros_Cilindros_Rotados = [];

	for j in 1:length(Arreglo_Centros_Cilindros)
		x = Mat_Cil_Z*Arreglo_Centros_Cilindros[j];
		push!(Arreglo_Centros_Cilindros_Rotados, x);
	end

	Arreglo_Posiciones = []; #Arreglo donde irán las coordenadas de puntos válidos.

	for i in 1:Iteraciones
		#Llaves para iniciar el código.
		Calcular = true;
		Intentos = 0;

		try
			while Calcular == true && Intentos < 1e4
				Intentos += 1;

				#En este punto coloquemos una partícula, sin importarnos si estamos dentro o fuera del cilindro.
				PrevX = Tipo(rand() - parse(Tipo, "0.5"));
				PrevY = Tipo(rand() - parse(Tipo, "0.5"));
				PrevZ = Tipo(rand() - parse(Tipo, "0.5"));
				Prev_Pos = [PrevX, PrevY, PrevZ];

				#Calculemos si ese punto es un punto válido.
				#Primero rotemos para que los cilindros queden verticales.
				Pos = Mat_Cil_Z*Prev_Pos;

				Llaves = Bool[];

				#Ahora coloquemos el origen en el centro de cada cilindro rotado y veamos si el punto cae dentro.
				for j in 1:length(Arreglo_Centros_Cilindros_Rotados)
					Posicion_Punto = Pos - Arreglo_Centros_Cilindros_Rotados[j]; #Punto de la partícula a prueba tras la rotación y traslación espacial del cilindro J al origen.
					((Posicion_Punto[3] < LadoCil/2 && Posicion_Punto[3] > -(LadoCil/2)) && ((Posicion_Punto[1])^2 + (Posicion_Punto[2])^2) < (RadioCil)^2) ? (push!(Llaves, true)) : (push!(Llaves, false));
				end

				any(Llaves) ? (push!(Arreglo_Posiciones, Prev_Pos); Calcular = false) : (nothing)

				Intentos == 1e4 ? throw("Demasiados Intentos sin encontrar un punto válido") : nothing;
			end
		catch
			nothing
		end
	end

	return Arreglo_Posiciones
end

function dibuja_Cilindros(EjeCil, Radio, LadoCaja, ix, iy, iz)

	Eje_Cil_Uni = EjeCil/norm(EjeCil); #Vector unitario y paralelo al eje del cilindro.
	Mat_Z_Cil = mat_Rot(EjeCil); #Matriz de rotación del eje Z al eje del cilindro.
	Mat_Cil_Z = (Mat_Z_Cil)^-1; #Matriz de rotación del eje del cilindro al eje Z.

	L = longitud_Cilindro(LadoCaja, Mat_Cil_Z); #Longitud del cilindro.
	θ = range(0.0, stop=2*π, length=10);
	ZZ = range(-L/2., stop=L/2., length=100);

	X = [];
	Y = [];
	Z = [];

	for j in θ
		for i in ZZ
			jj = (2*π)*rand();
			ii = rand(ZZ);
			xx = [Radio*cos(jj); Radio*sin(jj); ii];
			xx = Mat_Z_Cil*xx;
			xx = xx + [ix, iy, iz];

			if -LadoCaja/2. < xx[1] < LadoCaja/2. && -LadoCaja/2. < xx[2] < LadoCaja/2. && -LadoCaja/2. < xx[3] < LadoCaja/2.
				push!(X, xx[1]);
				push!(Y, xx[2]);
				push!(Z, xx[3]);
			end
		end
	end

	A = surface!(X,Y,Z, aspect_ratio = :equal, markersize = 1, key = false);

	return A
end

function interseccion_1(Radio1, Radio2, d)
	#Intersección entre dos círculos, uno en el origen y otro en [d,0]. Como hay dos intersecciones calcula los dos puntos.
	Tipo = typeof(d);
	Cero = parse(Tipo, "0");
	Dos = parse(Tipo, "2");
	
	X1 = (Radio1^2 - Radio2^2 + d^2)/(Dos*d);

	if (-Radio1^4 + Dos*Radio1^2*Radio2^2 - Radio2^4 + Dos*Radio1^2*d^2 + Dos*Radio2^2*d^2 - d^4) > Cero
		Y1 = -sqrt(-Radio1^4 + Dos*Radio1^2*Radio2^2 - Radio2^4 + Dos*Radio1^2*d^2 + Dos*Radio2^2*d^2 - d^4)/(Dos*d);
		return [X1; Y1], [X1; -Y1]
	else
		return [Inf; Inf], [Inf; Inf]
	end
end

function rotacion(X1, X2)
	#Matriz de rotación 2D para pasar de X1 a X2. También regresa la distancia d entre los puntos.
	d = norm(X1-X2);
	Tipo = typeof(d);
	Cero = parse(Tipo, "0");
	Uno = parse(Tipo, "1");

	E1 = [Uno; Cero];
	E2 = [Cero; Uno];

	CC = dot(X2 - X1, E1)/d;
	SS = dot(X2 - X1, E2)/d;

	return [[CC -SS]; [SS CC]], d
end

function interseccion_2(Radio1, Radio2, Centro1, Centro2)
	#Intersección entre 2 círculos de radios r1 y r2, y centros x1 y x2.

	#Solucionamos la intersección como si los círculos estuvieran en el origen y en el vector [d,0].
	R, d = rotacion(Centro1, Centro2);
	Prev_Colision1, Prev_Colision2 = interseccion_1(Radio1, Radio2, d);

	if Prev_Colision1[1] == Inf
		return [Inf; Inf], [Inf; Inf]
	end

	#Con la matriz de rotación regresamos las soluciones al caso en el que no están centrados los círculos en el origen y en el (0,d).
	Colision1 = R*Prev_Colision1 + Centro1;
	Colision2 = R*Prev_Colision2 + Centro1;

	return Colision1, Colision2
end

function tiempo_Trayectoria(Posicion_Inicial, Velocidad, Posicion_Final, Centro, Radio) ###No se usa la velocidad.
	#Tiempo que tarda en recorrerse de una posición inicial x, a una final xf.
	#La velocidad inicial es Velocidad, el centro de la circunferencia x1 y el radio Radio.

	Coseno = (Posicion_Inicial[1] - Centro[1])/Radio;
	Seno = (Posicion_Inicial[2] - Centro[2])/Radio;

	Tipo = typeof(Radio);
	Cero = parse(Tipo, "0");
	Uno = parse(Tipo, "1");
	Dos = parse(Tipo, "2");

	Cosωt = Coseno*(Posicion_Final[1] - Centro[1])/Radio + Seno*(Posicion_Final[2] - Centro[2])/Radio;
	Senωt = Coseno*(Posicion_Final[2] - Centro[2])/Radio - Seno*(Posicion_Final[1] - Centro[1])/Radio;

	if Cosωt > Uno
		Cosωt = Uno;

		return Cero
	end

	if Cosωt < -Uno
		Cosωt = -Uno;
	end

	if Senωt < Cero
		Tiempo = acos(Cosωt)*Radio;
		return Tiempo
	end

	if Senωt > Cero
		Tiempo = (Dos*π - acos(Cosωt))*Radio;
		return Tiempo
	end
end

function centro_Circunferencia(Velocidad, Posicion, Radio)
	#Centro de la circunferencia de radio R que da origen al movimiento circular con posición inicial Posicion y velocidad V.

	Velocidad = Velocidad/norm(Velocidad);
	Velocidad_Ortogonal = [Velocidad[2], -Velocidad[1]];
	Centro = Velocidad_Ortogonal*Radio + Posicion

	return Centro
end

function velocidad_Particula(Posicion, Centro)
	#Dada una posición y el centro de la circunferencia, calcula la velocidad de la partícula.

	Direccion = Posicion - Centro;
	Direccion = Direccion/norm(Direccion);

	Velocidad = [Direccion[2], -Direccion[1]]
end

function avance(Posicion_Inicial, Velocidad, Radio, Δt)
	#Avance de una partícula que se mueve sobre un círculo de radio R, a una velocidad inicial V, con posición inicial Pos_Ini y por un intervalo de tiempo Δt.

	Centro = centro_Circunferencia(Velocidad, Posicion_Inicial, Radio);
	Tipo = typeof(Posicion_Inicial[1]);
	ω = norm(Velocidad)/Radio;

	Posicion_Final = zeros(Tipo, 2);

	Coseno = (Posicion_Inicial[1] - Centro[1])/Radio;
	Seno = (Posicion_Inicial[2] - Centro[2])/Radio;

	Posicion_Final[1] = Radio*(cos(-ω*Δt)*Coseno - sin(-ω*Δt)*Seno) + Centro[1];
	Posicion_Final[2] = Radio*(cos(-ω*Δt)*Seno + sin(-ω*Δt)*Coseno) + Centro[2];

	Velocidad_Final = velocidad_Particula(Posicion_Final, Centro);

	return Posicion_Final, Velocidad_Final
end

function colision(Posicion_Inicial, Velocidad, Radio_Particula, Centro_Obstaculo, Radio_Obstaculo)
	#Calcula la posición de colisión y el tiempo que tarda para que esto suceda.
	Centro_Particula = centro_Circunferencia(Velocidad, Posicion_Inicial, Radio_Particula);

	Colision1, Colision2 = interseccion_2(Radio_Particula, Radio_Obstaculo, Centro_Particula, Centro_Obstaculo);

	if Colision2[1] == Inf || Colision2[1] == -Inf || Colision2[1] == NaN
		return [Inf; Inf], Inf
	end

	Tiempo = tiempo_Trayectoria(Posicion_Inicial, Velocidad, Colision2, Centro_Particula, Radio_Particula);
	return Colision2, Tiempo
end

function velocidad_Tras_Colision(Centro_Obstaculo, Posicion_Colision, Velocidad)
	#Velocidad de la partícula tras la colisión.

	Normal_Colision = (Centro_Obstaculo - Posicion_Colision);
	Tipo = typeof(Centro_Obstaculo[1]);

	Normal_Colision = Normal_Colision/norm(Normal_Colision);
	Velocidad_Normal_Colision = dot(Normal_Colision, Velocidad)*Normal_Colision;

	Velocidad = Velocidad - parse(Tipo, "2")*Velocidad_Normal_Colision;
	Velocidad = Velocidad/norm(Velocidad)

	return Velocidad
end

function colision_Condensada(Posicion_Inicial, Velocidad_Inicial, Radio_Particula, Centro_Obstaculo, Radio_Obstaculo)
	#Calcula el tiempo de colisión t, la posición de colisión y la velocidad final tras colision.

	Posicion_Colision, Tiempo_Colision = colision(Posicion_Inicial, Velocidad_Inicial, Radio_Particula, Centro_Obstaculo, Radio_Obstaculo);
	Centro_Particula = centro_Circunferencia(Velocidad_Inicial, Posicion_Inicial, Radio_Particula);
	Prev_Velocidad_Final = velocidad_Particula(Posicion_Colision, Centro_Particula);
	Velocidad_Final = velocidad_Tras_Colision(Centro_Obstaculo, Posicion_Colision, Prev_Velocidad_Final);

	return Tiempo_Colision, Posicion_Colision, Velocidad_Final
end


function colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, EjeCil, Centro_Cilindro, Radio_Cilindro, Longitud_Cilindro, Radio_Particula)

	Mat_Z_Cil = mat_Rot(EjeCil); #Matriz de rotación que transforma al eje Z en el eje del cilindro.
	Mat_Cil_Z = (Mat_Z_Cil)^-1; #Matriz de rotación que transforma al eje del cilindro en el eje Z.

	Tipo = typeof(Radio_Particula);
	Mat_3D_2D = [[parse(Tipo,"1") parse(Tipo,"0") parse(Tipo,"0")];[parse(Tipo,"0") parse(Tipo,"1") parse(Tipo,"0")]];

	Posicion_Cilindro_Vertical = Mat_Cil_Z*Posicion_Inicial;
	Velocidad_Cilindro_Vertical = Mat_Cil_Z*Velocidad_Inicial;
	Centro_Cilindro_Vertical = Mat_Cil_Z*Centro_Cilindro;
    h = Posicion_Cilindro_Vertical[3] - Centro_Cilindro_Vertical[3]
	#Veamos si el plano de la trayectoria de la partícula corta al cilindro, si no es el caso se regresa la posicion inicial, la velocidad inicial y un tiempo infinito.
	if abs(h) > Longitud_Cilindro/2
		return Posicion_Inicial, Velocidad_Inicial, Inf, 4, 4, h, [Inf, Inf]
	end
    
    
	Posicion_Cil_Vert_2D = Mat_3D_2D*Posicion_Cilindro_Vertical;
	Velocidad_Cil_Vert_2D = Mat_3D_2D*Velocidad_Cilindro_Vertical;
	Velocidad_Cil_Vert_2D = Velocidad_Cil_Vert_2D/norm(Velocidad_Cil_Vert_2D);
	Centro_Cil_Vert_2D = Mat_3D_2D*Centro_Cilindro_Vertical;

	Tiempo_Colision, Posicion_Colision, Velocidad_Colision = colision_Condensada(Posicion_Cil_Vert_2D, Velocidad_Cil_Vert_2D, Radio_Particula, Centro_Cil_Vert_2D, Radio_Cilindro);

	Posicion_Colision_Centrada_Cilindro = Posicion_Colision - Centro_Cil_Vert_2D;

	θ_Pos_Col_Cen_Cil = atan(Posicion_Colision_Centrada_Cilindro[2], Posicion_Colision_Centrada_Cilindro[1]);
	θ_Velocidad = atan(Velocidad_Colision[2], Velocidad_Colision[1]);
#    θ_Velocidad = θ_Pos_Col_Cen_Cil - θ_Velocidad
#	θ_Velocidad = mod(θ_Pos_Col_Cen_Cil - θ_Velocidad + π/2, π) - π/2;
#	θ_Pos_Col_Cen_Cil = mod(θ_Pos_Col_Cen_Cil, π);

	Posicion_Colision_3D = push!(Posicion_Colision, Posicion_Cilindro_Vertical[3]);
	Velocidad_Colision_3D = push!(Velocidad_Colision, Velocidad_Cilindro_Vertical[3]);
	Posicion_Colision_3D = Mat_Z_Cil*Posicion_Colision_3D;
	Velocidad_Colision_3D = Mat_Z_Cil*Velocidad_Colision_3D;

	Velocidad_Colision_3D = Velocidad_Colision_3D/norm(Velocidad_Colision_3D);

	return Posicion_Colision_3D, Velocidad_Colision_3D, Tiempo_Colision, θ_Pos_Col_Cen_Cil, θ_Velocidad, h, Centro_Cil_Vert_2D
end


function avanza(Posicion_Inicial, Velocidad_Inicial, EjeCil, Radio_Particula, Δt)
	Mat_Z_Cil = mat_Rot(EjeCil); #Matriz que transforma al eje Z en el eje del cilindro.
	Mat_Cil_Z = (Mat_Z_Cil)^-1; #Matriz que transforma al eje del cilindro en el eje Z.

	Tipo = typeof(Radio_Particula);
	Mat_3D_2D = [[parse(Tipo,"1") parse(Tipo,"0") parse(Tipo,"0")];[parse(Tipo,"0") parse(Tipo,"1") parse(Tipo,"0")]];

	Posicion_Cilindro_Vertical = Mat_Cil_Z*Posicion_Inicial;
	Velocidad_Cilindro_Vertical = Mat_Cil_Z*Velocidad_Inicial;

	Posicion_Cil_Vert_2D = Mat_3D_2D*Posicion_Cilindro_Vertical;
	Velocidad_Cil_Vert_2D = Mat_3D_2D*Velocidad_Cilindro_Vertical;
	Velocidad_Cil_Vert_2D = Velocidad_Cil_Vert_2D/norm(Velocidad_Cil_Vert_2D);

	Posicion_Intermedia, Velocidad_Intermedia = avance(Posicion_Cil_Vert_2D, Velocidad_Cil_Vert_2D, Radio_Particula, Δt);

	Posicion_Intermedia_3D = push!(Posicion_Intermedia, Posicion_Cilindro_Vertical[3]);
	Velocidad_Intermedia_3D = push!(Velocidad_Intermedia, Velocidad_Cilindro_Vertical[3]);

	Posicion_Intermedia_3D = Mat_Z_Cil*Posicion_Intermedia_3D;
	Velocidad_Intermedia_3D = Mat_Z_Cil*Velocidad_Intermedia_3D;
	Velocidad_Intermedia_3D = Velocidad_Intermedia_3D/norm(Velocidad_Intermedia_3D);

	return Posicion_Intermedia_3D, Velocidad_Intermedia_3D
end

function lorentz_Cuasi_Magnetico(Posicion_Inicial, Velocidad_Inicial, EjeCil, LadoCaja, Radio_Particula, Radio_Cilindro, Tiempo_Vuelo, Centros, N_Caja, Eficiencia = true)
	if Eficiencia
		Tipo = typeof(LadoCaja);
		EjeCil = EjeCil/norm(EjeCil)

		Mat_Z_Cil = mat_Rot(EjeCil);
		Mat_Cil_Z = (Mat_Z_Cil)^-1;

		L = longitud_Cilindro(LadoCaja, Mat_Cil_Z);
		N_Total = zeros(3);
		δT = parse(BigFloat, ".1");

		i = 0;
		Θc = [];
		Θv = [];
		Contador = 0;
		Tim = 0;

		while Tiempo_Vuelo > Tipo(0)
			i += 1;
			t = Inf;
			Posicion_Inicial, Coordenada_Caja = caja_Periodica(Posicion_Inicial, LadoCaja);
			N_Total += Coordenada_Caja;

			Posicion_1 = Posicion_Inicial;
			Velocidad_1 = Velocidad_Inicial;

			Centros_Cilindros = [Tipo(0), Tipo(0)];
			Coordenada_X_Cil = 0;
			Coordenada_Y_Cil = 0;

			for Coordenada_Z_Cil in -parse(Tipo, "1"):parse(Tipo, "1")
				Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, h, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, EjeCil, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, L, Radio_Particula);

				if Tiempo_Colision_1 < t #Calcula la colisión que consuma menos tiempo de vuelo (la primera colisión)
					t = Tiempo_Colision_1;
					Posicion_1 = copy(Posicion_2);
					Velocidad_1 = copy(Velocidad_2);
					Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);

					push!(Θc, θ_Pos_Col_Cen_Cil_1);
					push!(Θv, θ_Velocidad_1);

					Contador += 1;
				end
			end

			if t < δT
				Tiempo_Vuelo -= t;
				Posicion_Inicial = Posicion_1;
				Velocidad_Inicial = Velocidad_1;

				NN = Mat_Cil_Z*(N_Total + N_Caja);
				push!(Centros, Centros_Cilindros + [NN[1], NN[2]])
			else
				Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, EjeCil, Radio_Particula, δT)
				Tiempo_Vuelo -= δT;
				Tim += δT;
			end
		end

		return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, Centros
	else
		Tipo = typeof(LadoCaja);
		EjeCil = EjeCil/norm(EjeCil);

		Mat_Z_Cil = mat_Rot(EjeCil);
		Mat_Cil_Z = (Mat_Z_Cil)^-1;

		L = longitud_Cilindro(LadoCaja, Mat_Cil_Z);
		N_Total = zeros(3);
		δT = parse(BigFloat, ".1");

		i = 0;
		Θc = [];
		Θv = [];
		Contador = 0;
		Tim = 0;

		while Tiempo_Vuelo > Tipo(0)
			i += 1;
			t = Inf;
			Posicion_Inicial, Coordenada_Caja = caja_Periodica(Posicion_Inicial, LadoCaja);
			N_Total += Coordenada_Caja;

			Posicion_1 = Posicion_Inicial;
			Velocidad_1 = Velocidad_Inicial;

			Centros_Cilindros = [Tipo(0), Tipo(0)];

			for ix in -parse(Tipo, "1"):parse(Tipo, "1")
				for iy in -parse(Tipo, "1"):parse(Tipo, "1")
					for iz in -parse(Tipo, "1"):parse(Tipo, "1")

						Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, h, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, EjeCil, [ix, iy, iz], Radio_Cilindro, L, Radio_Particula);           
                        
						if Tiempo_Colision_1 < t

							t = Tiempo_Colision_1;
							Posicion_1 = copy(Posicion_2);
							Velocidad_1 = copy(Velocidad_2);
							Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);

							push!(Θc, θ_Pos_Col_Cen_Cil_1);
							push!(Θv, θ_Velocidad_1);

							Contador += 1;
						end

					end
				end
			end

			if t < δT
				Tiempo_Vuelo -= t;
				Posicion_Inicial = Posicion_1;
				Velocidad_Inicial = Velocidad_1;

				NN = Mat_Cil_Z*(N_Total + N_Caja);
				push!(Centros, Centros_Cilindros + [NN[1], NN[2]])
			else
				Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, EjeCil, Radio_Particula, δT)
				Tiempo_Vuelo -= δT;
				Tim += δT;
			end
		end

		return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, Centros

	end
end


function lorentz_Cuasi_Magnetico2(Posicion_Inicial, Velocidad_Inicial, EjeCil, LadoCaja, Radio_Particula, Radio_Cilindro, Tiempo_Vuelo, Centros, N_Caja, Eficiencia = true)
    δ = min(.5,Tiempo_Vuelo)
    δ2 = min(.1, Tiempo_Vuelo)
	if Eficiencia
		Tipo = typeof(LadoCaja);
		EjeCil = EjeCil/norm(EjeCil)

		Mat_Z_Cil = mat_Rot(EjeCil);
		Mat_Cil_Z = (Mat_Z_Cil)^-1;

		L = longitud_Cilindro(LadoCaja, Mat_Cil_Z);
		N_Total = zeros(3);
		δT = min(parse(BigFloat, "1000.1"), Tiempo_Vuelo);

		i = 0;
		Θc = [];
		Θv = [];
        H = [];
		Contador = 0;
#		Tim = 0;
        

		while Tiempo_Vuelo > Tipo(0)
			i += 1;
			t = Inf;
			Posicion_Inicial, Coordenada_Caja = caja_Periodica(Posicion_Inicial, LadoCaja);
			N_Total += Coordenada_Caja;
			Posicion_1 = Posicion_Inicial;
			Velocidad_1 = Velocidad_Inicial;
			Centros_Cilindros = [Tipo(0), Tipo(0)];
			Coordenada_X_Cil = parse(Tipo, "0");
			Coordenada_Y_Cil = parse(Tipo, "0");
            θ_Pos_Col_Cen_Cil_2 = 0;
            θ_Velocidad_2 = 0;
            h2 = 0
			for Coordenada_Z_Cil in [-parse(Tipo, "1"),parse(Tipo, "0"),parse(Tipo, "1")]
				Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, h, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, EjeCil, [Coordenada_X_Cil, Coordenada_Y_Cil, Coordenada_Z_Cil], Radio_Cilindro, L, Radio_Particula);
				if Tiempo_Colision_1 < t #Calcula la colisión que consuma menos tiempo de vuelo (la primera colisión)
                    nn = [BigFloat(0.),BigFloat(0.),Coordenada_Z_Cil]
					t = Tiempo_Colision_1;
					Posicion_1 = copy(Posicion_2);
					Velocidad_1 = copy(Velocidad_2);
					Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);
					Contador += 1;
                    θ_Velocidad_2 = copy(θ_Velocidad_1)
                    θ_Pos_Col_Cen_Cil_2 = copy(θ_Pos_Col_Cen_Cil_1)
                    h2 = h
				end
			end
                push!(H,h2)
                push!(Θc, θ_Pos_Col_Cen_Cil_2);
                push!(Θv, θ_Velocidad_2);

			if t < δT
				Tiempo_Vuelo -= t;
				Posicion_Inicial = Posicion_1;
				Velocidad_Inicial = Velocidad_1;
				NN = Mat_Cil_Z*(N_Total + N_Caja);
				push!(Centros, Centros_Cilindros + [NN[1], NN[2]])
			else
				Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, EjeCil, Radio_Particula, δ2)#
				Tiempo_Vuelo -= δ2;
#				Tim += 0.1;
			end
		end
		return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, H, Centros
	else
		Tipo = typeof(LadoCaja);
		EjeCil = EjeCil/norm(EjeCil);
		Mat_Z_Cil = mat_Rot(EjeCil);
		Mat_Cil_Z = (Mat_Z_Cil)^-1;
		L = longitud_Cilindro(LadoCaja, Mat_Cil_Z);
		N_Total = zeros(3);
		δT = min(parse(BigFloat, "1000.1"), Tiempo_Vuelo);
		i = 0;
		Θc = [];
		Θv = [];
        H = [];
		Contador = 0;
		Tim = 0;
		while Tiempo_Vuelo > Tipo(0)
            h2 = 0;
            θ_Pos_Col_Cen_Cil_2 = 0;
            θ_Velocidad_2 = 0;
			i += 1;
			t = Inf;
			Posicion_Inicial, Coordenada_Caja = caja_Periodica(Posicion_Inicial, LadoCaja);
			N_Total += Coordenada_Caja;

			Posicion_1 = Posicion_Inicial;
			Velocidad_1 = Velocidad_Inicial;
			Centros_Cilindros = [Tipo(0), Tipo(0)];
			for ix in -parse(Tipo, "1"):parse(Tipo, "1")
				for iy in -parse(Tipo, "1"):parse(Tipo, "1")
					for iz in -parse(Tipo, "1"):parse(Tipo, "1")
						Posicion_2, Velocidad_2, Tiempo_Colision_1, θ_Pos_Col_Cen_Cil_1, θ_Velocidad_1, h, Centro_Cilindro_Colision_2D = colision_Cilindro(Posicion_Inicial, Velocidad_Inicial, EjeCil, [ix, iy, iz], Radio_Cilindro, L, Radio_Particula);
						if Tiempo_Colision_1 < t
                            nn = [BigFloat(0.),BigFloat(0.),iz]
							t = Tiempo_Colision_1;
							Posicion_1 = copy(Posicion_2);
							Velocidad_1 = copy(Velocidad_2);
							Centros_Cilindros = copy(Centro_Cilindro_Colision_2D);
                            Contador += 1;
                            θ_Velocidad_2 = copy(θ_Velocidad_1)
                            θ_Pos_Col_Cen_Cil_2 = copy(θ_Pos_Col_Cen_Cil_1)
                            h2 = h
						end

					end
				end
			end
                push!(H,h2);
                push!(Θc, θ_Pos_Col_Cen_Cil_2);
                push!(Θv, θ_Velocidad_2);
			if t < δT
				Tiempo_Vuelo -= t;
				Posicion_Inicial = Posicion_1;
				Velocidad_Inicial = Velocidad_1;
				NN = Mat_Cil_Z*(N_Total + N_Caja);
				push!(Centros, Centros_Cilindros + [NN[1], NN[2]])
			else
                Posicion_Inicial, Velocidad_Inicial = avanza(Posicion_Inicial, Velocidad_Inicial, EjeCil, Radio_Particula, δ)
				Tiempo_Vuelo -= δ;
				Tim += δ;
			end
		end
		return Posicion_Inicial, Velocidad_Inicial, N_Total, Θc, Θv, H, Centros
	end
end



function dibuja_Circulo(Radio_Circulo, Centro_Circulo, color = "red")

	Tipo = typeof(Radio_Circulo);
	X = [Radio_Circulo];
	Y = [parse(Tipo, "0")];

	for i in 1:99
		x = Radio_Circulo*cos(2*π/100 * i);
		y = Radio_Circulo*sin(2*π/100 * i);

		push!(X, x);
		push!(Y, y);
	end

	push!(X, Radio_Circulo);
	push!(Y, parse(Tipo, "0"));

	X .+= Centro_Circulo[1];
	Y .+= Centro_Circulo[2];

	plot!(X, Y, color = color)
end

########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################

#= Θv = [];
Θc = [];

A = [];
Centros = [];
N_Particulas = 1000;
Tiempo_Vuelo_1 = 10000;
MSD = zeros(Tiempo_Vuelo_1);
LadoCaja = parse(BigFloat, "1.");
EjeCil = [BigFloat(1),BigFloat(2),π]/norm([BigFloat(1),BigFloat(2),π]);

Radio_Trayectoria = parse(BigFloat, "7.");
Radio_Cilindro = BigFloat(.15);
Mat_Z_Cil = mat_Rot(EjeCil);
Mat_Cil_Z = (Mat_Z_Cil)^-1;

for j in 1:N_Particulas
	tic()
	println(j)

	velocidad_inicial_1 = vel_Rand(Mat_Z_Cil);
	velocidad_inicial_1 /= norm(velocidad_inicial_1);

	pos_inicial = posiciones_Iniciales(EjeCil, Radio_Cilindro, LadoCaja, 1);
	pos_inicial_1 = pos_inicial[1];

	X_pos_inicial_1 = [pos_inicial_1[1]];
	Y_pos_inicial_1 = [pos_inicial_1[2]];
	Z_pos_inicial_1 = [pos_inicial_1[3]];

	pos_inicial_1_cilindro_vertical = Mat_Cil_Z*pos_inicial_1;

	X_Cil_Vert = [pos_inicial_1_cilindro_vertical[1]];
	Y_Cil_Vert = [pos_inicial_1_cilindro_vertical[2]];
	Z_Cil_Vert = [pos_inicial_1_cilindro_vertical[3]];

	Time = [0];
	Tiempo_De_Vuelo = BigFloat(1);
	N_Total = zeros(BigFloat, 3);
	Centros = [];

	for i in 1:Tiempo_Vuelo_1
		pos_inicial_1, velocidad_inicial_1, n, Θc, Θv, Centros = lorentz_Cuasi_Magnetico(pos_inicial_1, velocidad_inicial_1, EjeCil, LadoCaja, Radio_Trayectoria, Radio_Cilindro, Tiempo_De_Vuelo, Centros, N_Total)
		N_Total += n;

		pos_inicial_1_cilindro_vertical = Mat_Cil_Z*(pos_inicial_1 + N_Total);

		push!(X_Cil_Vert, pos_inicial_1_cilindro_vertical[1]);
		push!(Y_Cil_Vert, pos_inicial_1_cilindro_vertical[2]);

		MSD[i] += X_Cil_Vert[i]^2 + Y_Cil_Vert[i]^2;
	end

	toc()
end

MSD /= N_Particulas
=#
import math

def dividir_circunferencia(r):
    # loop para dividir a circunferência em 12 fatias iguais
    for i in range(1,360,30):
        # ângulo da fatia em relação ao eixo x
        angulo = i * math.pi/180
        
        # coordenadas do ponto na circunferência no formato do xml
        print(f"<point x=\"{round(r * math.cos(angulo),3)}\" y=\"{round(r * math.sin(angulo),3)}\" z=\"0\"/>")   

dividir_circunferencia(21)




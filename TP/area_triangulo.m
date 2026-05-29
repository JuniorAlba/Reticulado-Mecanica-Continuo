% Area con signo de un triangulo (positiva = orientacion antihoraria)
function A_sig = area_triangulo(xi, yi, xj, yj, xk, yk)
    A_sig = 0.5 * ((xj-xi)*(yk-yi) - (xk-xi)*(yj-yi));
end

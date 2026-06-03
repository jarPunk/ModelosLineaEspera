"""
Logica de calculo para modelos de lineas de espera.
"""

from __future__ import annotations
import math
import re

SIGNIFICADOS: dict[str, str] = {
    "P0": "probabilidad de que no haya unidades en el sistema",
    "Pk = Pb": "probabilidad de bloqueo (todas las k posiciones ocupadas)",
    "Pw": "probabilidad de que una llegada tenga que esperar",
    "ρ": "utilizacion del sistema o de los servidores",
    "a = λ/μ": "trafico ofrecido",
    "Lq": "numero promedio de unidades en cola",
    "L": "numero promedio de unidades en el sistema",
    "Wq": "tiempo promedio en cola",
    "W": "tiempo promedio en el sistema",
    "E[S]": "tiempo medio de servicio",
    "E[S^2]": "segundo momento del tiempo de servicio",
    "λe": "tasa de llegada efectiva",
    "ρ = λe/(k·μ)": "utilizacion promedio por servidor",
}

def _etiqueta_con_significado(simbolo: str) -> str:
    if simbolo in SIGNIFICADOS:
        return f"{simbolo} ({SIGNIFICADOS[simbolo]})"

    match_pn = re.fullmatch(r"P(\d+)", simbolo)
    if match_pn:
        return f"{simbolo} (probabilidad de tener n={match_pn.group(1)} unidades en el sistema)"

    return simbolo

def _a_flotante(texto: str) -> float:
    texto_limpio = texto.strip().replace(",", ".")
    if not texto_limpio:
        raise ValueError("Ingresa un valor numerico.")

    if "/" in texto_limpio:
        partes = texto_limpio.split("/")
        if len(partes) != 2:
            raise ValueError("Formato de fraccion invalido.")

        numerador = float(partes[0].strip())
        denominador = float(partes[1].strip())
        if abs(denominador) < 1e-12:
            raise ValueError("El denominador no puede ser 0.")
        return numerador / denominador

    return float(texto_limpio)

def _asegurar_mayor(nombre: str, valor: float, limite: float = 0.0) -> None:
    if valor <= limite:
        raise ValueError(f"{nombre} debe ser mayor que {limite}.")

def _asegurar_entero(nombre: str, valor: float, minimo: int = 1) -> int:
    entero = int(valor)
    if abs(valor - entero) > 1e-9:
        raise ValueError(f"{nombre} debe ser entero.")
    if entero < minimo:
        raise ValueError(f"{nombre} debe ser >= {minimo}.")
    return entero

def calc_mm1(lambd: float, mu: float, n_opcional: int | None = None) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    rho = lambd / mu
    if rho >= 1: raise ValueError("Sistema inestable: λ < μ.")
    p0 = 1 - rho
    lq = (rho**2) / (1 - rho)
    l = rho / (1 - rho)
    wq = lq / lambd
    w = l / lambd
    resultados = {"P0": p0, "ρ": rho, "Lq": lq, "L": l, "Wq": wq, "W": w}
    if n_opcional is not None: resultados[f"P{n_opcional}"] = (rho**n_opcional) * p0
    return resultados

def calc_md1(lambd: float, mu: float) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    rho = lambd / mu
    if rho >= 1: raise ValueError("Sistema inestable: λ < μ.")
    p0 = 1 - rho
    pw = rho
    lq = (rho**2) / (2 * (1 - rho))
    wq = lq / lambd
    w = wq + (1 / mu)
    l = lambd * w
    return {"P0": p0, "Pw": pw, "ρ": rho, "Lq": lq, "L": l, "Wq": wq, "W": w}

def _pn_mmk(a: float, p0: float, k: int, n: int) -> float:
    if n < 0: raise ValueError("n debe ser >= 0.")
    if n < k: return ((a**n) / math.factorial(n)) * p0
    return ((a**n) / (math.factorial(k) * (k ** (n - k)))) * p0

def calc_mmk(lambd: float, mu: float, k: int, n_opcional: int | None = None) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if k < 1: raise ValueError("k debe ser >= 1.")
    a = lambd / mu
    rho = lambd / (k * mu)
    if rho >= 1: raise ValueError("Sistema inestable: λ < k·μ.")
    suma = sum((a**n) / math.factorial(n) for n in range(k))
    termino_cola = (a**k) / (math.factorial(k) * (1 - rho))
    p0 = 1 / (suma + termino_cola)
    p_espera = termino_cola * p0
    lq = p0 * ((a**k) * rho) / (math.factorial(k) * ((1 - rho) ** 2))
    wq = lq / lambd
    w = wq + (1 / mu)
    l = lambd * w
    resultados = {"P0": p0, "Pw": p_espera, "ρ": rho, "a = λ/μ": a, "Lq": lq, "L": l, "Wq": wq, "W": w}
    if n_opcional is not None: resultados[f"P{n_opcional}"] = _pn_mmk(a, p0, k, n_opcional)
    return resultados

def calc_mg1(lambd: float, mu: float, cs2: float) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if cs2 < 0: raise ValueError("Cs^2 debe ser >= 0.")
    es = 1 / mu
    es2 = (1 + cs2) / (mu**2)
    rho = lambd * es
    if rho >= 1: raise ValueError("Sistema inestable: λ·E[S] < 1.")
    wq = (lambd * es2) / (2 * (1 - rho))
    w = wq + es
    lq = lambd * wq
    l = lambd * w
    p0 = 1 - rho
    return {"P0": p0, "ρ": rho, "E[S]": es, "E[S^2]": es2, "Lq": lq, "L": l, "Wq": wq, "W": w}

def _erlang_b(a: float, k: int) -> float:
    b = 1.0
    for i in range(1, k + 1):
        b = (a * b) / (i + (a * b))
    return b

def calc_mgk_bloqueo(lambd: float, mu: float, k: int, j_opcional: int | None = None) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if k < 1: raise ValueError("k debe ser >= 1.")
    a = lambd / mu
    pb = _erlang_b(a, k)
    p0 = 1 / sum((a**n) / math.factorial(n) for n in range(k + 1))
    lambda_efectiva = lambd * (1 - pb)
    l = a * (1 - pb)
    utilizacion_media = l / k
    w = 1 / mu
    resultados = {"P0": p0, "Pk = Pb": pb, "Pw": a, "λe": lambda_efectiva, "ρ": utilizacion_media, "Lq": 0.0, "L": l, "Wq": 0.0, "W": w}
    if j_opcional is not None: resultados[f"P{j_opcional}"] = ((a**j_opcional) / math.factorial(j_opcional)) * p0
    return resultados

def calc_mm1_fuente_finita(lambd: float, mu: float, n_fuentes: int, n_opcional: int | None = None) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if n_fuentes < 1: raise ValueError("N debe ser >= 1.")
    razon = lambd / mu
    suma = sum((math.factorial(n_fuentes) / math.factorial(n_fuentes - n)) * (razon**n) for n in range(n_fuentes + 1))
    p0 = 1 / suma
    lq = n_fuentes - ((lambd + mu) / lambd) * (1 - p0)
    l = lq + (1 - p0)
    lambda_efectiva = (n_fuentes - l) * lambd
    if lambda_efectiva <= 0: raise ValueError("Tasa efectiva de llegada <= 0.")
    wq = lq / lambda_efectiva
    w = wq + (1 / mu)
    pw = 1 - p0
    resultados = {"P0": p0, "Pw": pw, "Lq": lq, "L": l, "Wq": wq, "W": w}
    if n_opcional is not None: resultados[f"P{n_opcional}"] = (math.factorial(n_fuentes) / math.factorial(n_fuentes - n_opcional)) * (razon**n_opcional) * p0
    return resultados

def calc_mmk_infinito_n(lambd: float, mu: float, k: int, n_fuentes: int) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if k < 1: raise ValueError("k debe ser >= 1.")
    if n_fuentes < 1: raise ValueError("N debe ser >= 1.")
    coeficientes = [1.0]
    for n in range(1, n_fuentes + 1):
        tasa_llegada_prev = (n_fuentes - (n - 1)) * lambd
        tasa_servicio_n = min(n, k) * mu
        coeficientes.append(coeficientes[-1] * (tasa_llegada_prev / tasa_servicio_n))
    suma = sum(coeficientes)
    p0 = 1 / suma
    p = [p0 * c for c in coeficientes]
    l = sum(n * p[n] for n in range(n_fuentes + 1))
    lambda_efectiva = sum((n_fuentes - n) * lambd * p[n] for n in range(n_fuentes + 1))
    w = l / lambda_efectiva if lambda_efectiva > 0 else 0.0
    rho = lambda_efectiva / (k * mu)
    wq = max(0.0, w - (1 / mu))
    lq = lambda_efectiva * wq
    limite = min(k - 1, n_fuentes)
    p_espera = 1 - sum(p[n] for n in range(limite + 1))
    return {"P0": p0, "Pw": p_espera, "λe": lambda_efectiva, "ρ = λe/(k·μ)": rho, "Lq": lq, "L": l, "Wq": wq, "W": w}

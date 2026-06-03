"""
Calculadora grafica de modelos de lineas de espera.

Modelos incluidos:
- M/M/1
- M/D/1
- M/M/k (Erlang C)
- M/G/1 (Pollaczek-Khinchine)
- M/G/K con clientes bloqueados eliminados (Erlang B)
- M/M/1 con fuente finita
- M/M/k/infinito/N con fuente finita
"""

from __future__ import annotations

import math
import re
import tkinter as tk
from tkinter import ttk, messagebox


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

SALIDAS_POR_MODELO: dict[str, list[str]] = {
    "M/M/1": ["P0", "ρ", "Lq", "L", "Wq", "W", "Pn (opcional)"],
    "M/D/1": ["P0", "Pw", "ρ", "Lq", "L", "Wq", "W"],
    "M/M/k": ["P0", "Pw", "ρ", "a = λ/μ", "Lq", "L", "Wq", "W", "Pn (opcional)"],
    "M/G/1": ["P0", "ρ", "E[S]", "E[S^2]", "Lq", "L", "Wq", "W"],
    "M/G/K": ["P0", "Pk = Pb", "a = λ/μ", "λe", "ρ", "Lq", "L", "Wq", "W", "Pj (opcional)"],
    "M/M/1 Fuente Finita": ["P0", "Pw", "Lq", "L", "Wq", "W", "Pn (opcional)"],
    "M/M/k/∞/N": ["P0", "Pw", "λe", "ρ = λe/(k·μ)", "Lq", "L", "Wq", "W"],
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
            raise ValueError("Formato de fraccion invalido. Usa: numerador/denominador.")

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
    if n_opcional is not None and n_opcional < 0:
        raise ValueError("n debe ser >= 0.")

    rho = lambd / mu
    if rho >= 1:
        raise ValueError("Sistema inestable: se requiere λ < μ.")

    p0 = 1 - rho
    lq = (rho**2) / (1 - rho)
    l = rho / (1 - rho)
    wq = lq / lambd
    w = l / lambd

    resultados = {
        "P0": p0,
        "ρ": rho,
        "Lq": lq,
        "L": l,
        "Wq": wq,
        "W": w,
    }

    if n_opcional is not None:
        resultados[f"P{n_opcional}"] = (rho**n_opcional) * p0

    return resultados


def calc_md1(lambd: float, mu: float) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)

    rho = lambd / mu
    if rho >= 1:
        raise ValueError("Sistema inestable: se requiere λ < μ.")

    p0 = 1 - rho
    pw = rho
    lq = (rho**2) / (2 * (1 - rho))
    wq = lq / lambd
    w = wq + (1 / mu)
    l = lambd * w

    return {
        "P0": p0,
        "Pw": pw,
        "ρ": rho,
        "Lq": lq,
        "L": l,
        "Wq": wq,
        "W": w,
    }


def _pn_mmk(a: float, p0: float, k: int, n: int) -> float:
    if n < 0:
        raise ValueError("n debe ser >= 0.")
    if n < k:
        return ((a**n) / math.factorial(n)) * p0
    return ((a**n) / (math.factorial(k) * (k ** (n - k)))) * p0


def calc_mmk(lambd: float, mu: float, k: int, n_opcional: int | None = None) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if k < 1:
        raise ValueError("k debe ser >= 1.")

    a = lambd / mu
    rho = lambd / (k * mu)
    if rho >= 1:
        raise ValueError("Sistema inestable: se requiere λ < k·μ.")

    suma = sum((a**n) / math.factorial(n) for n in range(k))
    termino_cola = (a**k) / (math.factorial(k) * (1 - rho))
    p0 = 1 / (suma + termino_cola)

    p_espera = termino_cola * p0
    lq = p0 * ((a**k) * rho) / (math.factorial(k) * ((1 - rho) ** 2))
    wq = lq / lambd
    w = wq + (1 / mu)
    l = lambd * w

    resultados = {
        "P0": p0,
        "Pw": p_espera,
        "ρ": rho,
        "a = λ/μ": a,
        "Lq": lq,
        "L": l,
        "Wq": wq,
        "W": w,
    }

    if n_opcional is not None:
        resultados[f"P{n_opcional}"] = _pn_mmk(a, p0, k, n_opcional)

    return resultados


def calc_mg1(lambd: float, mu: float, cs2: float) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if cs2 < 0:
        raise ValueError("Cs^2 debe ser >= 0.")

    es = 1 / mu
    es2 = (1 + cs2) / (mu**2)
    rho = lambd * es

    if rho >= 1:
        raise ValueError("Sistema inestable: se requiere λ·E[S] < 1.")

    wq = (lambd * es2) / (2 * (1 - rho))
    w = wq + es
    lq = lambd * wq
    l = lambd * w
    p0 = 1 - rho

    return {
        "P0": p0,
        "ρ": rho,
        "E[S]": es,
        "E[S^2]": es2,
        "Lq": lq,
        "L": l,
        "Wq": wq,
        "W": w,
    }


def _erlang_b(a: float, k: int) -> float:
    b = 1.0
    for i in range(1, k + 1):
        b = (a * b) / (i + (a * b))
    return b


def calc_mgk_bloqueo(lambd: float, mu: float, k: int, j_opcional: int | None = None) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if k < 1:
        raise ValueError("k debe ser >= 1.")
    if j_opcional is not None and (j_opcional < 0 or j_opcional > k):
        raise ValueError("j debe estar entre 0 y k.")

    a = lambd / mu
    pb = _erlang_b(a, k)
    p0 = 1 / sum((a**n) / math.factorial(n) for n in range(k + 1))
    lambda_efectiva = lambd * (1 - pb)
    l = a * (1 - pb)
    utilizacion_media = l / k
    w = 1 / mu

    resultados = {
        "P0": p0,
        "Pk = Pb": pb,
        "a = λ/μ": a,
        "λe": lambda_efectiva,
        "ρ": utilizacion_media,
        "Lq": 0.0,
        "L": l,
        "Wq": 0.0,
        "W": w,
    }

    if j_opcional is not None:
        resultados[f"P{j_opcional}"] = ((a**j_opcional) / math.factorial(j_opcional)) * p0

    return resultados


def calc_mm1_fuente_finita(lambd: float, mu: float, n_fuentes: int, n_opcional: int | None = None) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if n_fuentes < 1:
        raise ValueError("N debe ser >= 1.")
    if n_opcional is not None and (n_opcional < 0 or n_opcional > n_fuentes):
        raise ValueError("n debe estar entre 0 y N.")

    razon = lambd / mu

    suma = 0.0
    for n in range(n_fuentes + 1):
        suma += (math.factorial(n_fuentes) / math.factorial(n_fuentes - n)) * (razon**n)
    p0 = 1 / suma

    lq = n_fuentes - ((lambd + mu) / lambd) * (1 - p0)
    l = lq + (1 - p0)

    lambda_efectiva = (n_fuentes - l) * lambd
    if lambda_efectiva <= 0:
        raise ValueError("No se puede calcular Wq: la tasa efectiva de llegada es <= 0.")

    wq = lq / lambda_efectiva
    w = wq + (1 / mu)
    pw = 1 - p0

    resultados = {
        "P0": p0,
        "Pw": pw,
        "Lq": lq,
        "L": l,
        "Wq": wq,
        "W": w,
    }

    if n_opcional is not None:
        pn = (math.factorial(n_fuentes) / math.factorial(n_fuentes - n_opcional)) * (razon**n_opcional) * p0
        resultados[f"P{n_opcional}"] = pn

    return resultados


def calc_mmk_infinito_n(lambd: float, mu: float, k: int, n_fuentes: int) -> dict[str, float]:
    _asegurar_mayor("λ", lambd)
    _asegurar_mayor("μ", mu)
    if k < 1:
        raise ValueError("k debe ser >= 1.")
    if n_fuentes < 1:
        raise ValueError("N debe ser >= 1.")

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

    return {
        "P0": p0,
        "Pw": p_espera,
        "λe": lambda_efectiva,
        "ρ = λe/(k·μ)": rho,
        "Lq": lq,
        "L": l,
        "Wq": wq,
        "W": w,
    }


class QueueApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Calculadora de Colas")
        self.geometry("1020x700")
        self.minsize(940, 620)

        self.configure(bg="#EEF3F8")
        self._configurar_estilos()
        self._construir_ui()

    def _configurar_estilos(self) -> None:
        style = ttk.Style(self)
        style.theme_use("clam")

        style.configure("TNotebook", background="#EEF3F8", borderwidth=0)
        style.configure(
            "TNotebook.Tab",
            font=("Segoe UI", 10, "bold"),
            padding=(14, 10),
            background="#D9E3EF",
            foreground="#243447",
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", "#2E5B8A")],
            foreground=[("selected", "#FFFFFF")],
        )

        style.configure("Card.TFrame", background="#FFFFFF")
        style.configure("Panel.TFrame", background="#F7FAFD")
        style.configure("Card.TLabel", background="#FFFFFF", foreground="#243447", font=("Segoe UI", 10))
        style.configure("Panel.TLabel", background="#F7FAFD", foreground="#324A5E", font=("Segoe UI", 10))
        style.configure("PanelTitle.TLabel", background="#F7FAFD", foreground="#1F3B5C", font=("Segoe UI", 10, "bold"))
        style.configure("Title.TLabel", background="#FFFFFF", foreground="#163A59", font=("Segoe UI", 20, "bold"))
        style.configure("Subtitle.TLabel", background="#FFFFFF", foreground="#54708A", font=("Segoe UI", 10))
        style.configure("Run.TButton", font=("Segoe UI", 10, "bold"), padding=(14, 9), background="#2E5B8A", foreground="#FFFFFF")
        style.map("Run.TButton", background=[("active", "#244A72")])
        style.configure("TEntry", fieldbackground="#FFFFFF", bordercolor="#BDCAD8", lightcolor="#BDCAD8", darkcolor="#BDCAD8")

    def _construir_ui(self) -> None:
        header = ttk.Frame(self, style="Card.TFrame", padding=16)
        header.pack(fill="x", padx=18, pady=(16, 10))

        ttk.Label(header, text="Modelos de Lineas de Espera", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            header,
            text="Ingresa variables en cada pestaña y presiona Calcular. Notacion: P0, Pn/Pk, Lq, L, Wq, W, ρ, λe.",
            style="Subtitle.TLabel",
        ).pack(anchor="w", pady=(4, 0))

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=18, pady=(0, 18))

        self._crear_tab_modelo(
            titulo="M/M/1",
            campos=[
                ("lambd", "λ Tasa de llegada"),
                ("mu", "μ Tasa de servicio"),
                ("n", "n para calcular Pn (opcional)"),
            ],
            calculo=self._calc_tab_mm1,
        )

        self._crear_tab_modelo(
            titulo="M/D/1",
            campos=[
                ("lambd", "λ Tasa media de llegadas"),
                ("mu", "μ Tasa media de servicios"),
            ],
            calculo=self._calc_tab_md1,
        )

        self._crear_tab_modelo(
            titulo="M/M/k",
            campos=[
                ("lambd", "λ Tasa de llegada"),
                ("mu", "μ Tasa de servicio"),
                ("k", "k Servidores"),
                ("n", "n para calcular Pn (opcional)"),
            ],
            calculo=self._calc_tab_mmk,
        )

        self._crear_tab_modelo(
            titulo="M/G/1",
            campos=[("lambd", "λ Tasa de llegada"), ("mu", "μ Tasa media de servicio"), ("cs2", "Cs^2 Variabilidad del servicio")],
            calculo=self._calc_tab_mg1,
        )

        self._crear_tab_modelo(
            titulo="M/G/K",
            campos=[
                ("lambd", "λ Tasa de llegada"),
                ("mu", "μ Tasa media de servicio"),
                ("k", "k Canales"),
                ("j", "j para calcular Pj (opcional)"),
            ],
            calculo=self._calc_tab_mgkk,
        )

        self._crear_tab_modelo(
            titulo="M/M/1 Fuente Finita",
            campos=[
                ("lambd", "λ Tasa de llegada por unidad"),
                ("mu", "μ Tasa de servicio"),
                ("n", "N Tamano de la poblacion"),
                ("n_prob", "n para calcular Pn (opcional)"),
            ],
            calculo=self._calc_tab_mm1_fuente_finita,
        )

        self._crear_tab_modelo(
            titulo="M/M/k/∞/N",
            campos=[
                ("lambd", "λ Tasa por fuente"),
                ("mu", "μ Tasa de servicio"),
                ("k", "k Servidores"),
                ("n", "N Poblacion finita"),
            ],
            calculo=self._calc_tab_mmk_infinito_n,
        )

    def _crear_tab_modelo(self, titulo: str, campos: list[tuple[str, str]], calculo) -> None:
        tab = ttk.Frame(self.notebook, style="Card.TFrame", padding=18)
        self.notebook.add(tab, text=titulo)

        ttk.Label(tab, text=f"Configuracion del modelo {titulo}", style="Subtitle.TLabel").pack(anchor="w", pady=(0, 10))

        frm_superior = ttk.Frame(tab, style="Card.TFrame")
        frm_superior.pack(fill="x")

        frm_inputs = ttk.Frame(frm_superior, style="Panel.TFrame", padding=12)
        frm_inputs.pack(side="left", fill="x", expand=True)

        frm_info = ttk.Frame(frm_superior, style="Panel.TFrame", padding=12)
        frm_info.pack(side="right", fill="y", padx=(18, 0))

        salidas = SALIDAS_POR_MODELO.get(titulo, [])
        texto_salidas = "\n".join(f"- {_etiqueta_con_significado(s)}" for s in salidas) if salidas else "- Sin datos"
        ttk.Label(frm_info, text="Resultados disponibles", style="PanelTitle.TLabel").pack(anchor="w")
        ttk.Label(
            frm_info,
            text=texto_salidas,
            style="Panel.TLabel",
            justify="left",
            wraplength=330,
        ).pack(anchor="w", pady=(4, 0))

        entries: dict[str, ttk.Entry] = {}
        for fila, (key, etiqueta) in enumerate(campos):
            ttk.Label(frm_inputs, text=etiqueta, style="Panel.TLabel").grid(row=fila, column=0, sticky="w", pady=6, padx=(0, 12))
            entrada = ttk.Entry(frm_inputs, width=28)
            entrada.grid(row=fila, column=1, sticky="w", pady=6)
            entries[key] = entrada

        frm_costos = ttk.Frame(frm_superior, style="Panel.TFrame", padding=12)
        frm_costos.pack(side="right", fill="y", padx=(18, 0))

        ttk.Label(frm_costos, text="Costo total (opcional)", style="PanelTitle.TLabel").grid(row=0, column=0, columnspan=2, sticky="w")
        ttk.Label(frm_costos, text="C_unidad (costo por unidad)", style="Panel.TLabel").grid(row=1, column=0, sticky="w", pady=(8, 4), padx=(0, 10))
        entry_c_unidad = ttk.Entry(frm_costos, width=16)
        entry_c_unidad.grid(row=1, column=1, sticky="w", pady=(8, 4))

        ttk.Label(frm_costos, text="C_servidor (costo por servidor)", style="Panel.TLabel").grid(row=2, column=0, sticky="w", pady=4, padx=(0, 10))
        entry_c_servidor = ttk.Entry(frm_costos, width=16)
        entry_c_servidor.grid(row=2, column=1, sticky="w", pady=4)

        ttk.Label(frm_costos, text="Metrica base", style="Panel.TLabel").grid(row=3, column=0, sticky="w", pady=4, padx=(0, 10))
        metrica_var = tk.StringVar(value="L")
        combo_metrica = ttk.Combobox(frm_costos, textvariable=metrica_var, values=("L", "Lq"), width=13, state="readonly")
        combo_metrica.grid(row=3, column=1, sticky="w", pady=4)

        ttk.Label(
            frm_costos,
            text="CT = C_unidad * metrica + C_servidor * k",
            style="Panel.TLabel",
            wraplength=250,
            justify="left",
        ).grid(row=4, column=0, columnspan=2, sticky="w", pady=(8, 0))

        btn = ttk.Button(tab, text="Calcular", style="Run.TButton")
        btn.pack(anchor="w", pady=(14, 10))

        cont_salida = ttk.Frame(tab, style="Panel.TFrame", padding=8)
        cont_salida.pack(fill="both", expand=True)

        salida = tk.Text(
            cont_salida,
            height=14,
            bg="#FDFEFF",
            fg="#203040",
            insertbackground="#203040",
            font=("Consolas", 10),
            relief="solid",
            bd=1,
            padx=10,
            pady=10,
        )
        scroll = ttk.Scrollbar(cont_salida, orient="vertical", command=salida.yview)
        salida.configure(yscrollcommand=scroll.set)
        salida.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")

        def on_click() -> None:
            try:
                resultado = calculo(entries)
                salida.config(state="normal")
                salida.delete("1.0", tk.END)
                for nombre, valor in resultado.items():
                    etiqueta = _etiqueta_con_significado(nombre)
                    salida.insert(tk.END, f"{etiqueta:<78} = {valor:.6f}\n")

                c_unidad_txt = entry_c_unidad.get().strip()
                c_servidor_txt = entry_c_servidor.get().strip()
                if c_unidad_txt or c_servidor_txt:
                    if not c_unidad_txt or not c_servidor_txt:
                        raise ValueError("Para calcular CT, ingresa ambos costos: C_unidad y C_servidor.")

                    c_unidad = _a_flotante(c_unidad_txt)
                    c_servidor = _a_flotante(c_servidor_txt)
                    if c_unidad < 0 or c_servidor < 0:
                        raise ValueError("Los costos no pueden ser negativos.")

                    metrica_ct = metrica_var.get().strip()
                    if metrica_ct not in resultado:
                        raise ValueError(f"La metrica {metrica_ct} no esta disponible para este modelo.")

                    if "k" in entries and entries["k"].get().strip():
                        s = _asegurar_entero("k", _a_flotante(entries["k"].get()), minimo=1)
                    else:
                        s = 1

                    ct = (c_unidad * float(resultado[metrica_ct])) + (c_servidor * s)
                    salida.insert(tk.END, "\n")
                    salida.insert(tk.END, "-" * 98 + "\n")
                    salida.insert(tk.END, "Costo Total (CT)\n")
                    salida.insert(tk.END, f"CT = ({c_unidad:.6f} * {metrica_ct}) + ({c_servidor:.6f} * k={s})\n")
                    salida.insert(tk.END, f"CT = {ct:.6f}\n")

                salida.config(state="disabled")
            except Exception as exc:
                messagebox.showerror("Dato invalido", str(exc))

        btn.configure(command=on_click)

    def _leer_entrada(self, entries: dict[str, ttk.Entry], key: str) -> float:
        valor = _a_flotante(entries[key].get())
        return valor

    def _calc_tab_mm1(self, entries: dict[str, ttk.Entry]) -> dict[str, float]:
        lambd = self._leer_entrada(entries, "lambd")
        mu = self._leer_entrada(entries, "mu")
        n_texto = entries["n"].get().strip()
        if n_texto:
            try:
                n_opcional = int(n_texto)
            except ValueError as exc:
                raise ValueError("n debe ser un entero >= 0.") from exc
            if n_opcional < 0:
                raise ValueError("n debe ser >= 0.")
        else:
            n_opcional = None
        return calc_mm1(lambd, mu, n_opcional=n_opcional)

    def _calc_tab_md1(self, entries: dict[str, ttk.Entry]) -> dict[str, float]:
        lambd = self._leer_entrada(entries, "lambd")
        mu = self._leer_entrada(entries, "mu")
        return calc_md1(lambd, mu)

    def _calc_tab_mmk(self, entries: dict[str, ttk.Entry]) -> dict[str, float]:
        lambd = self._leer_entrada(entries, "lambd")
        mu = self._leer_entrada(entries, "mu")
        k = _asegurar_entero("k", self._leer_entrada(entries, "k"), minimo=1)
        n_texto = entries["n"].get().strip()
        if n_texto:
            try:
                n_opcional = int(n_texto)
            except ValueError as exc:
                raise ValueError("n debe ser un entero >= 0.") from exc
            if n_opcional < 0:
                raise ValueError("n debe ser >= 0.")
        else:
            n_opcional = None
        return calc_mmk(lambd, mu, k, n_opcional=n_opcional)

    def _calc_tab_mg1(self, entries: dict[str, ttk.Entry]) -> dict[str, float]:
        lambd = self._leer_entrada(entries, "lambd")
        mu = self._leer_entrada(entries, "mu")
        cs2 = self._leer_entrada(entries, "cs2")
        return calc_mg1(lambd, mu, cs2)

    def _calc_tab_mgkk(self, entries: dict[str, ttk.Entry]) -> dict[str, float]:
        lambd = self._leer_entrada(entries, "lambd")
        mu = self._leer_entrada(entries, "mu")
        k = _asegurar_entero("k", self._leer_entrada(entries, "k"), minimo=1)
        j_texto = entries["j"].get().strip()
        if j_texto:
            try:
                j_opcional = int(j_texto)
            except ValueError as exc:
                raise ValueError("j debe ser un entero entre 0 y k.") from exc
        else:
            j_opcional = None
        return calc_mgk_bloqueo(lambd, mu, k, j_opcional=j_opcional)

    def _calc_tab_mm1_fuente_finita(self, entries: dict[str, ttk.Entry]) -> dict[str, float]:
        lambd = self._leer_entrada(entries, "lambd")
        mu = self._leer_entrada(entries, "mu")
        n_fuentes = _asegurar_entero("N", self._leer_entrada(entries, "n"), minimo=1)

        n_texto = entries["n_prob"].get().strip()
        if n_texto:
            try:
                n_opcional = int(n_texto)
            except ValueError as exc:
                raise ValueError("n debe ser un entero entre 0 y N.") from exc
        else:
            n_opcional = None

        return calc_mm1_fuente_finita(lambd, mu, n_fuentes, n_opcional=n_opcional)

    def _calc_tab_mmk_infinito_n(self, entries: dict[str, ttk.Entry]) -> dict[str, float]:
        lambd = self._leer_entrada(entries, "lambd")
        mu = self._leer_entrada(entries, "mu")
        k = _asegurar_entero("k", self._leer_entrada(entries, "k"), minimo=1)
        n = _asegurar_entero("N", self._leer_entrada(entries, "n"), minimo=1)
        return calc_mmk_infinito_n(lambd, mu, k, n)


def ejecutar() -> None:
    app = QueueApp()
    app.mainloop()


if __name__ == "__main__":
    ejecutar()

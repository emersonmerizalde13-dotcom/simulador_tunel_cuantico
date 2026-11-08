
"use strict";

/* ============================================
   SIMULADOR CUÁNTICO UNEMI
   Modelo estacionario + paquete de ondas 1D
============================================ */

// CONSTANTES FÍSICAS
const MASA = 9.1093837e-31;
const HBAR = 1.054571817e-34;
const EV = 1.602176634e-19;
const NM = 1e-9;

// Unidades del paquete: nm, fs y eV
const HBAR_EV_FS = 0.6582119569;
const CINETICA_EV_NM2 = 0.03809982;

// CONFIGURACIÓN INICIAL
const INICIAL = {
  E: 1,
  V0: 2,
  a: 0.5
};

// CONFIGURACIÓN DEL PAQUETE
const N = 1024;
const XMIN = -20;
const XMAX = 20;
const L = XMAX - XMIN;
const DX = L / N;

const DT_NUMERICO = 0.02; // fs
// La simulación no tiene un límite fijo de tiempo.

// ESTADO
let modelo = null;
let paquete = null;
let fase = 0;
let reproduciendo = false;
let ultimoTiempo = null;
let acumulador = 0;

const $ = (id) => document.getElementById(id);

const controles = {
  energia: $("energia"),
  barrera: $("barrera"),
  ancho: $("ancho"),
  real: $("mostrarReal"),
  imag: $("mostrarImag"),
  densidad: $("mostrarDensidad")
};

const graficas = {
  energia: $("grafEnergia"),
  onda: $("grafOnda"),
  densidad: $("grafDensidad")
};

/* ============================================
   1. NAVEGACIÓN
============================================ */

function mostrarPantalla(id) {
  ["bienvenida", "creditos", "simulador"].forEach(
    (nombre) => {
      $(nombre).classList.toggle(
        "oculto",
        nombre !== id
      );
    }
  );

  if (id !== "simulador") {
    reproduciendo = false;
    ultimoTiempo = null;
  }

  window.scrollTo(0, 0);

  if (id === "simulador") {
    requestAnimationFrame(dibujarTodo);
  }
}

$("empezar").addEventListener("click", () => {
  mostrarPantalla("creditos");
});

$("volverInicio").addEventListener("click", () => {
  mostrarPantalla("bienvenida");
});

$("entrar").addEventListener("click", () => {
  mostrarPantalla("simulador");
});

$("regresar").addEventListener("click", () => {
  mostrarPantalla("bienvenida");
});

/* ============================================
   2. NÚMEROS COMPLEJOS
============================================ */

function complejo(re, im = 0) {
  return { re, im };
}

function sumar(z, w) {
  return complejo(z.re + w.re, z.im + w.im);
}

function restar(z, w) {
  return complejo(z.re - w.re, z.im - w.im);
}

function multiplicar(z, w) {
  return complejo(
    z.re * w.re - z.im * w.im,
    z.re * w.im + z.im * w.re
  );
}

function dividir(z, w) {
  const d = w.re * w.re + w.im * w.im;

  return complejo(
    (z.re * w.re + z.im * w.im) / d,
    (z.im * w.re - z.re * w.im) / d
  );
}

function modulo2(z) {
  return z.re * z.re + z.im * z.im;
}

function expI(angulo) {
  return complejo(
    Math.cos(angulo),
    Math.sin(angulo)
  );
}

/* ============================================
   3. MODELO ESTACIONARIO
============================================ */

function calcularModelo(E, V0, a) {
  const energiaJ = E * EV;
  const barreraJ = V0 * EV;
  const anchoM = a * NM;

  const k = Math.sqrt(
    2 * MASA * energiaJ
  ) / HBAR;

  const diferencia = energiaJ - barreraJ;

  let tipo;
  let parametro;
  let D;

  if (diferencia < -1e-14 * EV) {
    tipo = "tunel";

    const kappa = Math.sqrt(
      -2 * MASA * diferencia
    ) / HBAR;

    parametro = kappa;

    const u = kappa * anchoM;
    const B = (k * k - kappa * kappa) /
      (2 * k * kappa);

    D = complejo(
      Math.cosh(u),
      -B * Math.sinh(u)
    );

  } else if (diferencia > 1e-14 * EV) {
    tipo = "superior";

    const q = Math.sqrt(
      2 * MASA * diferencia
    ) / HBAR;

    parametro = q;

    const u = q * anchoM;
    const B = (k * k + q * q) /
      (2 * k * q);

    D = complejo(
      Math.cos(u),
      -B * Math.sin(u)
    );

  } else {
    tipo = "igual";
    parametro = 0;

    D = complejo(
      1,
      -k * anchoM / 2
    );
  }

  const t = dividir(
    expI(-k * anchoM),
    D
  );

  const T = Math.max(
    0,
    Math.min(1, modulo2(t))
  );

  const R = 1 - T;

  const psiA = multiplicar(
    t,
    expI(k * anchoM)
  );

  function ondaInterior(xM) {
    const d = xM - anchoM;
    let factor;

    if (tipo === "tunel") {
      const u = parametro * d;

      factor = complejo(
        Math.cosh(u),
        (k / parametro) * Math.sinh(u)
      );

    } else if (tipo === "superior") {
      const u = parametro * d;

      factor = complejo(
        Math.cos(u),
        (k / parametro) * Math.sin(u)
      );

    } else {
      factor = complejo(1, k * d);
    }

    return multiplicar(psiA, factor);
  }

  const r = restar(
    ondaInterior(0),
    complejo(1)
  );

  function onda(xNm) {
    const xM = xNm * NM;

    if (xM < 0) {
      return sumar(
        expI(k * xM),
        multiplicar(
          r,
          expI(-k * xM)
        )
      );
    }

    if (xM <= anchoM) {
      return ondaInterior(xM);
    }

    return multiplicar(
      t,
      expI(k * xM)
    );
  }

  return {
    E, V0, a, T, R, tipo, onda
  };
}

/* ============================================
   4. TRANSFORMADA RÁPIDA DE FOURIER
============================================ */

function fft(re, im, inversa = false) {
  const n = re.length;

  // Reordenamiento de bits
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;

    while (j & bit) {
      j ^= bit;
      bit >>= 1;
    }

    j ^= bit;

    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }

  for (let len = 2; len <= n; len <<= 1) {
    const angulo =
      (inversa ? 2 : -2) * Math.PI / len;

    const wr = Math.cos(angulo);
    const wi = Math.sin(angulo);

    for (let i = 0; i < n; i += len) {
      let ur = 1;
      let ui = 0;

      for (let j = 0; j < len / 2; j++) {
        const p = i + j;
        const q = p + len / 2;

        const vr = re[q] * ur - im[q] * ui;
        const vi = re[q] * ui + im[q] * ur;

        const ar = re[p];
        const ai = im[p];

        re[p] = ar + vr;
        im[p] = ai + vi;

        re[q] = ar - vr;
        im[q] = ai - vi;

        const nr = ur * wr - ui * wi;
        const ni = ur * wi + ui * wr;

        ur = nr;
        ui = ni;
      }
    }
  }

  if (inversa) {
    for (let i = 0; i < n; i++) {
      re[i] /= n;
      im[i] /= n;
    }
  }
}

/* ============================================
   5. CREAR PAQUETE DE ONDAS
============================================ */

function crearPaquete() {
  const E = Number(controles.energia.value);
  const V0 = Number(controles.barrera.value);
  const a = Number(controles.ancho.value);

  const re = new Float64Array(N);
  const im = new Float64Array(N);

  // Posición inicial y anchura espacial
  const x0 = -8;
  const sigma = 1.5;

  // Número de onda central
  const k0 = Math.sqrt(
    E / CINETICA_EV_NM2
  );

  let norma = 0;

  for (let i = 0; i < N; i++) {
    const x = XMIN + i * DX;

    const envolvente = Math.exp(
      -((x - x0) ** 2) /
      (4 * sigma * sigma)
    );

    const angulo = k0 * (x - x0);

    re[i] = envolvente * Math.cos(angulo);
    im[i] = envolvente * Math.sin(angulo);

    norma += (
      re[i] * re[i] +
      im[i] * im[i]
    ) * DX;
  }

  const factor = 1 / Math.sqrt(norma);

  for (let i = 0; i < N; i++) {
    re[i] *= factor;
    im[i] *= factor;
  }

  paquete = {
    re,
    im,
    E,
    V0,
    a,
    tiempo: 0
  };

  acumulador = 0;
}

/* ============================================
   6. EVOLUCIÓN DE SCHRÖDINGER
   Método split-step Fourier
============================================ */

function evolucionarPaquete(dt) {
  if (!paquete) return;

  const { re, im, V0, a } = paquete;

  // Medio paso de potencial
  aplicarPotencial(re, im, V0, a, dt / 2);

  // Pasar al espacio de momentos
  fft(re, im, false);

  // Paso cinético completo
  for (let i = 0; i < N; i++) {
    const indice = i < N / 2 ? i : i - N;

    const k = 2 * Math.PI * indice / L;

    const energiaK =
      CINETICA_EV_NM2 * k * k;

    const angulo =
      -energiaK * dt / HBAR_EV_FS;

    const c = Math.cos(angulo);
    const s = Math.sin(angulo);

    const r = re[i];
    const v = im[i];

    re[i] = r * c - v * s;
    im[i] = r * s + v * c;
  }

  // Volver al espacio de posiciones
  fft(re, im, true);

  // Segundo medio paso de potencial
  aplicarPotencial(re, im, V0, a, dt / 2);

  paquete.tiempo += dt;
}

function aplicarPotencial(re, im, V0, a, dt) {
  for (let i = 0; i < N; i++) {
    const x = XMIN + i * DX;

    const V =
      x >= 0 && x <= a ? V0 : 0;

    const angulo =
      -V * dt / HBAR_EV_FS;

    const c = Math.cos(angulo);
    const s = Math.sin(angulo);

    const r = re[i];
    const v = im[i];

    re[i] = r * c - v * s;
    im[i] = r * s + v * c;
  }
}

/* ============================================
   7. PREPARAR GRÁFICAS
============================================ */

function prepararCanvas(canvas) {
  const rect = canvas.getBoundingClientRect();

  const w = Math.max(300, rect.width);
  const h = Math.max(150, rect.height);

  const escala = Math.min(
    window.devicePixelRatio || 1,
    2
  );

  canvas.width = Math.round(w * escala);
  canvas.height = Math.round(h * escala);

  const ctx = canvas.getContext("2d");

  ctx.setTransform(
    escala, 0, 0, escala, 0, 0
  );

  ctx.clearRect(0, 0, w, h);

  return { ctx, w, h };
}

function crearPlano(
  canvas, xmin, xmax, ymin, ymax
) {
  const { ctx, w, h } = prepararCanvas(canvas);

  const izquierda = 52;
  const derecha = w - 18;
  const superior = 18;
  const inferior = h - 35;

  const px = (x) =>
    izquierda +
    (x - xmin) / (xmax - xmin) *
    (derecha - izquierda);

  const py = (y) =>
    inferior -
    (y - ymin) / (ymax - ymin) *
    (inferior - superior);

  ctx.save();

  ctx.strokeStyle = "rgba(150,190,220,0.15)";
  ctx.lineWidth = 1;

  for (let i = 0; i <= 5; i++) {
    const x = izquierda +
      (derecha - izquierda) * i / 5;

    ctx.beginPath();
    ctx.moveTo(x, superior);
    ctx.lineTo(x, inferior);
    ctx.stroke();
  }

  for (let i = 0; i <= 4; i++) {
    const y = superior +
      (inferior - superior) * i / 4;

    ctx.beginPath();
    ctx.moveTo(izquierda, y);
    ctx.lineTo(derecha, y);
    ctx.stroke();
  }

  ctx.strokeStyle = "#6385a8";
  ctx.beginPath();
  ctx.moveTo(izquierda, superior);
  ctx.lineTo(izquierda, inferior);
  ctx.lineTo(derecha, inferior);
  ctx.stroke();

  ctx.fillStyle = "#a9c2d9";
  ctx.font = "11px Arial";

  ctx.textAlign = "center";

  for (let i = 0; i <= 5; i++) {
    const x = xmin + (xmax - xmin) * i / 5;

    ctx.fillText(
      x.toFixed(1),
      px(x),
      inferior + 17
    );
  }

  ctx.textAlign = "right";

  for (let i = 0; i <= 4; i++) {
    const y = ymin + (ymax - ymin) * i / 4;

    ctx.fillText(
      y.toFixed(1),
      izquierda - 7,
      py(y) + 4
    );
  }

  ctx.fillText("x (nm)", derecha, h - 5);

  ctx.restore();

  return {
    ctx, px, py,
    xmin, xmax, ymin, ymax,
    izquierda, derecha,
    superior, inferior
  };
}

function dibujarBarrera(plano, a) {
  const { ctx, px, superior, inferior } = plano;

  ctx.save();

  ctx.fillStyle = "rgba(167,139,250,0.17)";
  ctx.fillRect(
    px(0),
    superior,
    px(a) - px(0),
    inferior - superior
  );

  ctx.strokeStyle = "#b48cff";
  ctx.setLineDash([5, 4]);

  [0, a].forEach((x) => {
    ctx.beginPath();
    ctx.moveTo(px(x), superior);
    ctx.lineTo(px(x), inferior);
    ctx.stroke();
  });

  ctx.restore();
}

function trazarCurva(
  plano, funcion, color, grosor = 2
) {
  const { ctx, px, py, xmin, xmax } = plano;

  ctx.save();
  ctx.beginPath();

  for (let i = 0; i <= 600; i++) {
    const x = xmin + (xmax - xmin) * i / 600;
    const y = funcion(x);

    if (!Number.isFinite(y)) continue;

    if (i === 0) {
      ctx.moveTo(px(x), py(y));
    } else {
      ctx.lineTo(px(x), py(y));
    }
  }

  ctx.strokeStyle = color;
  ctx.lineWidth = grosor;
  ctx.stroke();
  ctx.restore();
}

/* ============================================
   8. GRÁFICA 1: ENERGÍA POTENCIAL
============================================ */

function dibujarEnergia() {
  const { E, V0, a } = modelo;

  const ymax = Math.max(E, V0) * 1.3;

  const plano = crearPlano(
    graficas.energia,
    -2, a + 2, 0, ymax
  );

  dibujarBarrera(plano, a);

  const { ctx, px, py } = plano;

  ctx.save();

  // Barrera rectangular
  ctx.strokeStyle = "#b48cff";
  ctx.lineWidth = 2.5;

  ctx.beginPath();
  ctx.moveTo(px(-2), py(0));
  ctx.lineTo(px(0), py(0));
  ctx.lineTo(px(0), py(V0));
  ctx.lineTo(px(a), py(V0));
  ctx.lineTo(px(a), py(0));
  ctx.lineTo(px(a + 2), py(0));
  ctx.stroke();

  // Energía del electrón
  ctx.strokeStyle = "#00d9f5";
  ctx.setLineDash([7, 5]);

  ctx.beginPath();
  ctx.moveTo(px(-2), py(E));
  ctx.lineTo(px(a + 2), py(E));
  ctx.stroke();

  ctx.setLineDash([]);
  ctx.font = "12px Arial";

  ctx.fillStyle = "#00d9f5";
  ctx.textAlign = "left";

  ctx.fillText(
    "E = " + E.toFixed(2) + " eV",
    px(-1.9),
    Math.max(30, py(E) - 8)
  );

  ctx.fillStyle = "#c8adff";

  ctx.fillText(
    "V₀ = " + V0.toFixed(2) + " eV",
    px(0) + 6,
    Math.max(16, py(V0) - 8)
  );

  ctx.restore();
}

/* ============================================
   9. GRÁFICA 2: FUNCIÓN DE ONDA
============================================ */

function dibujarOnda() {
  const xmin = -2;
  const xmax = modelo.a + 2;

  let maximo = 1;

  for (let i = 0; i <= 300; i++) {
    const x = xmin + (xmax - xmin) * i / 300;

    maximo = Math.max(
      maximo,
      Math.sqrt(modulo2(modelo.onda(x)))
    );
  }

  const limite = maximo * 1.2;

  const plano = crearPlano(
    graficas.onda,
    xmin, xmax, -limite, limite
  );

  dibujarBarrera(plano, modelo.a);

  const rotacion = expI(-fase);

  if (controles.real.checked) {
    trazarCurva(
      plano,
      (x) => multiplicar(
        modelo.onda(x),
        rotacion
      ).re,
      "#00d9f5",
      2.2
    );
  }

  if (controles.imag.checked) {
    trazarCurva(
      plano,
      (x) => multiplicar(
        modelo.onda(x),
        rotacion
      ).im,
      "#ffb65c",
      2
    );
  }
}

/* ============================================
   10. GRÁFICA 3: PAQUETE DINÁMICO
============================================ */

function dibujarDensidad() {
  if (!paquete) return;

  const plano = crearPlano(
    graficas.densidad,
    XMIN, XMAX,
    0, 0.55
  );

  dibujarBarrera(plano, paquete.a);

  const { ctx, px, py } = plano;

  if (controles.densidad.checked) {
    ctx.save();

    // Relleno bajo la densidad
    ctx.beginPath();
    ctx.moveTo(px(XMIN), py(0));

    for (let i = 0; i < N; i++) {
      const x = XMIN + i * DX;

      const densidad =
        paquete.re[i] ** 2 +
        paquete.im[i] ** 2;

      ctx.lineTo(px(x), py(densidad));
    }

    ctx.lineTo(px(XMAX), py(0));
    ctx.closePath();

    ctx.fillStyle = "rgba(40,230,184,0.15)";
    ctx.fill();

    // Curva de densidad
    ctx.beginPath();

    for (let i = 0; i < N; i++) {
      const x = XMIN + i * DX;

      const densidad =
        paquete.re[i] ** 2 +
        paquete.im[i] ** 2;

      if (i === 0) {
        ctx.moveTo(px(x), py(densidad));
      } else {
        ctx.lineTo(px(x), py(densidad));
      }
    }

    ctx.strokeStyle = "#28e6b8";
    ctx.lineWidth = 2.5;
    ctx.stroke();

    ctx.restore();
  }

  // Tiempo de simulación
  ctx.save();

  ctx.fillStyle = "#d1e9f7";
  ctx.font = "12px Arial";
  ctx.textAlign = "right";

  ctx.fillText(
    "Tiempo: " +
    paquete.tiempo.toFixed(2) +
    " fs",
    plano.derecha,
    plano.superior + 14
  );

  ctx.restore();
}

function dibujarTodo() {
  if (!modelo || !paquete) return;

  dibujarEnergia();
  dibujarOnda();
  dibujarDensidad();
}

/* ============================================
   11. ACTUALIZAR RESULTADOS
============================================ */

function actualizar() {
  const E = Number(controles.energia.value);
  const V0 = Number(controles.barrera.value);
  const a = Number(controles.ancho.value);

  modelo = calcularModelo(E, V0, a);

  $("energiaValor").textContent =
    E.toFixed(2) + " eV";

  $("barreraValor").textContent =
    V0.toFixed(2) + " eV";

  $("anchoValor").textContent =
    a.toFixed(2) + " nm";

  $("transmision").textContent =
    (modelo.T * 100).toFixed(6) + " %";

  $("reflexion").textContent =
    (modelo.R * 100).toFixed(6) + " %";

  $("barraT").style.width =
    (modelo.T * 100) + "%";

  $("barraR").style.width =
    (modelo.R * 100) + "%";

  if (modelo.tipo === "tunel") {
    $("regimen").textContent =
      "E < V₀ · Túnel cuántico";

    $("explicacion").textContent =
      "La energía del electrón es menor " +
      "que la barrera. Aun así, existe " +
      "una probabilidad de transmisión " +
      "debido al efecto túnel cuántico.";

  } else if (modelo.tipo === "superior") {
    $("regimen").textContent =
      "E > V₀ · Dispersión cuántica";

    $("explicacion").textContent =
      "El electrón tiene energía superior " +
      "a la barrera, pero todavía puede " +
      "reflejarse por interferencia cuántica.";

  } else {
    $("regimen").textContent =
      "E = V₀ · Caso límite";

    $("explicacion").textContent =
      "La energía coincide con la altura " +
      "de la barrera. Se aplica la " +
      "solución límite del modelo.";
  }

  dibujarTodo();
}

/* ============================================
   12. CONTROLES
============================================ */

[
  controles.energia,
  controles.barrera,
  controles.ancho
].forEach((control) => {
  control.addEventListener("input", () => {
    reproduciendo = false;
    ultimoTiempo = null;
    fase = 0;

    crearPaquete();
    actualizar();
  });
});

[
  controles.real,
  controles.imag,
  controles.densidad
].forEach((control) => {
  control.addEventListener(
    "change",
    dibujarTodo
  );
});


$("iniciar").addEventListener("click", () => {
  reproduciendo = true;
  ultimoTiempo = null;
});


$("pausar").addEventListener("click", () => {
  reproduciendo = false;
  ultimoTiempo = null;
});

$("reiniciar").addEventListener("click", () => {
  reproduciendo = false;
  ultimoTiempo = null;
  fase = 0;

  crearPaquete();
  dibujarTodo();
});

$("reiniciarTodo").addEventListener("click", () => {
  reproduciendo = false;
  ultimoTiempo = null;
  fase = 0;

  controles.energia.value = INICIAL.E;
  controles.barrera.value = INICIAL.V0;
  controles.ancho.value = INICIAL.a;

  controles.real.checked = true;
  controles.imag.checked = true;
  controles.densidad.checked = true;

  crearPaquete();
  actualizar();
});

/* ============================================
   13. CICLO DE ANIMACIÓN
============================================ */


function animar(tiempo) {
  if (reproduciendo) {
    if (ultimoTiempo !== null) {
      // Segundos reales transcurridos entre fotogramas
      const dtReal = Math.min(
        (tiempo - ultimoTiempo) / 1000,
        0.05
      );

      // Animación de la segunda gráfica
      fase = (
        fase + dtReal * 2.5
      ) % (2 * Math.PI);

      // 1 segundo real = 2.5 femtosegundos simulados
      const velocidadSimulacion = 2.5;

      acumulador += dtReal * velocidadSimulacion;

      let pasos = 0;

      // Evolución física del paquete de ondas
      while (
        acumulador >= DT_NUMERICO &&
        pasos < 30
      ) {
        evolucionarPaquete(DT_NUMERICO);

        acumulador -= DT_NUMERICO;
        pasos++;
      }

      // Actualizar las gráficas
      if (
        !$("simulador").classList.contains("oculto")
      ) {
        dibujarOnda();
        dibujarDensidad();
      }
    }

    ultimoTiempo = tiempo;
  } else {
    ultimoTiempo = null;
  }

  requestAnimationFrame(animar);
}


/* ============================================
   14. INICIALIZACIÓN
============================================ */

window.addEventListener("resize", () => {
  if (
    !$("simulador").classList.contains("oculto")
  ) {
    dibujarTodo();
  }
});

// Crear paquete inicial
crearPaquete();

// Calcular y mostrar resultados
actualizar();

// Mostrar bienvenida al iniciar
mostrarPantalla("bienvenida");

// Activar el motor de animación
requestAnimationFrame(animar);

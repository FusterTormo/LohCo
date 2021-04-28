/*Cambiar la pestana que se va a mostrar, dependiendo del boton que se ha clicado*/
function abrir(pestana, boton) {
  var i, tabs, buts;
  tabs = document.getElementsByTagName("section");
  buts = document.getElementsByClassName("boton");
  for (i = 0; i < tabs.length; i++)
    tabs[i].className = "invisible";
  for (i = 0; i < buts.length; i++)
    buts[i].className = "boton"
  document.getElementById(pestana).className = "visible";
  boton.className = "boton select";
}

/*Buscar en la tabla de variantes mientras se escribe*/
function buscar(que) {
  var crom, pos, gen, tabla, filas, it, celda, cnt;
  tabla = document.getElementById("tabVar");
  filas = tabla.getElementsByTagName("tr");
  if (que == "posicion") {
    //Buscar por cromosoma-posicion
    crom = "chr" + document.getElementById("chr").value.toUpperCase();
    pos = document.getElementById("position").value;
    //Buscar solo por la posicion si el cromosoma esta vacio
    if (crom == "chr") {
      for (it = 0; it < filas.length; it ++) {
        celda = filas[it].getElementsByTagName("td")[1];
        if (celda) {
          cnt = celda.textContent || celda.innerText;
          if (cnt == pos)
            filas[it].style.display = "";
          else
            filas[it].style.display = "none";
        }
      }
    }
    // Buscar solo por el cromosoma si la posicion esta vacia
    if (pos == "") {
      for (it = 0; it < filas.length; it ++) {
        celda = filas[it].getElementsByTagName("td")[0];
        if (celda) {
          cnt = celda.textContent || celda.innerText;
          if (cnt == crom)
            filas[it].style.display = "";
          else
            filas[it].style.display = "none";
        }
      }
    }
    // Buscar por cromosoma y posicion
    if (crom != "chr" && pos != "") {
      for (it = 0; it < filas.length; it ++) {
        celda = filas[it].getElementsByTagName("td")[0];
        celda2 = filas[it].getElementsByTagName("td")[1];
        if (celda && celda2) {
          cnt = celda.textContent || celda.innerText;
          cnt2 = celda2.textContent || celda2.innerText;
          if (cnt == crom && cnt2 == pos)
            filas[it].style.display = "";
          else
            filas[it].style.display = "none";
        }
      }
    }
    // Mostrar todas las filas de la tabla en caso de que ambos campos esten vacios
    if (crom == "chr" && pos == "") {
      for (it = 0; it < filas.length; it ++) {
        filas[it].style.display = "";
      }
    }
  }
  else if (que == "gen") {
    // Buscar por gen
    gen = document.getElementById("gene").value.toUpperCase();
    for (it = 0; it < filas.length; it++) {
      celda = filas[it].getElementsByTagName("td")[4];
      if (celda) {
        cnt = celda.textContent || celda.innerText;
        if (cnt.indexOf(gen) > -1)
          filas[it].style.display = "";
        else
          filas[it].style.display = "none";
      }
    }
  }
}

/*Crea una ventana donde se muestra la imagen a la que se ha clicado para verla en grande*/
function crearModal() {
	var modal = document.getElementById("elModal");
	var imgs = document.getElementsByTagName("img");
	var contenido = document.getElementById("imagen");
	var it;
	for (it = 0; it < imgs.length; it++) {
		imgs[it].onclick = function() {
			modal.style.display = "block";
			contenido.src = this.src;
		}
	}

	var cerrar = document.getElementsByClassName("close")[0];
	cerrar.onclick = function() {
		modal.style.display = "none";
	}
}

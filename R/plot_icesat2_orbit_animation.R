#' Plot ICESat-2 orbital tracks as a 3D globe animation
#'
#' @description
#' Creates an interactive HTML animation of ICESat-2 orbital tracks from one or
#' more KML files or from a \code{terra::SpatVector} returned by
#' \code{rgt_extract()}.
#'
#' @param kml_files Character vector. Path to one or more KML files. If
#'   \code{NULL} and \code{rgt = NULL}, all example KML files available in
#'   \code{inst/extdata} are used.
#' @param rgt A \code{terra::SpatVector} returned by \code{rgt_extract()}.
#'   If provided, this object is used instead of \code{kml_files}.
#' @param output_dir Character. Folder where the HTML animation and required
#'   assets will be written. Default is \code{tempdir()}.
#' @param output_file Character. Name of the output HTML file.
#' @param track_speed Numeric. Initial animation speed. Default is \code{2}.
#' @param earth_rotation_speed Numeric. Initial Earth rotation speed. Default is \code{2}.
#' @param launch Logical. If \code{TRUE}, starts a local server and opens the HTML.
#'
#' @return Invisibly returns the path to the generated HTML file.
#'
#' @examples
#' \dontrun{
#' plot_icesat2_orbit_animation()
#'
#' rgt_line <- rgt_extract(h5 = atl03_h5, line = TRUE)
#'
#' plot_icesat2_orbit_animation(
#'   rgt = rgt_line,
#'   launch = TRUE
#' )
#' }
#'
#' @export
plot_icesat2_orbit_animation <- function(
    kml_files = NULL,
    rgt = NULL,
    output_dir = tempdir(),
    output_file = "ICESat2_orbit_animation.html",
    track_speed = 2,
    earth_rotation_speed = 2,
    launch = TRUE
) {

  requireNamespace("sf")
  requireNamespace("dplyr")
  requireNamespace("stringr")
  requireNamespace("jsonlite")
  requireNamespace("servr")
  requireNamespace("terra")

  earth_texture <- "Stylized_World_Topo_5400x2700.jpeg"

  extdata_dir <- system.file(
    "extdata",
    package = "ICESat2VegR"
  )

  if (extdata_dir == "") {
    stop("Could not find the package extdata folder.")
  }

  earth_texture_path <- file.path(
    extdata_dir,
    earth_texture
  )

  if (!file.exists(earth_texture_path)) {
    stop("Earth texture not found in inst/extdata: ", earth_texture)
  }

  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }

  file.copy(
    earth_texture_path,
    file.path(output_dir, earth_texture),
    overwrite = TRUE
  )

  all_tracks <- list()

  if (!is.null(rgt)) {

    if (!inherits(rgt, "SpatVector")) {
      stop("rgt must be a terra::SpatVector returned by rgt_extract().")
    }

    rgt_sf <- sf::st_as_sf(rgt)
    coords <- sf::st_coordinates(rgt_sf)

    track_df <- data.frame(
      track_id = "rgt_extract",
      track_order = 1L,
      step = seq_len(nrow(coords)),
      time_label = NA_character_,
      lon = coords[, "X"],
      lat = coords[, "Y"]
    )

  } else {

    if (is.null(kml_files)) {
      kml_files <- list.files(
        extdata_dir,
        pattern = "\\.kml$",
        full.names = TRUE
      )
    }

    if (length(kml_files) == 0) {
      stop("No KML files provided, no rgt object supplied, and no example KML files found in inst/extdata.")
    }

    if (any(!file.exists(kml_files))) {
      stop("Some KML files do not exist.")
    }

    kml_files <- sort(kml_files)

    for (i in seq_along(kml_files)) {

      kf <- kml_files[i]
      layers <- sf::st_layers(kf)$name

      layer_tracks <- list()

      for (layer_name in layers) {

        x <- sf::st_read(
          kf,
          layer = layer_name,
          quiet = TRUE
        )

        coords <- sf::st_coordinates(x)

        one_layer <- data.frame(
          track_id = tools::file_path_sans_ext(basename(kf)),
          track_order = i,
          step = seq_len(nrow(coords)),
          time_label = stringr::str_extract(layer_name, "\\d{2}:\\d{2}:\\d{2}"),
          lon = coords[, "X"],
          lat = coords[, "Y"]
        )

        layer_tracks[[layer_name]] <- one_layer
      }

      one_track <- dplyr::bind_rows(layer_tracks)

      if (nrow(one_track) > 2) {

        first_pt <- one_track[1, ]
        last_pt <- one_track[nrow(one_track), ]

        close_dist <- sqrt(
          (first_pt$lon - last_pt$lon)^2 +
            (first_pt$lat - last_pt$lat)^2
        )

        if (close_dist < 0.01) {
          one_track <- one_track[-nrow(one_track), ]
        }
      }

      one_track <- one_track |>
        dplyr::mutate(
          step = dplyr::row_number()
        )

      all_tracks[[i]] <- one_track
    }

    track_df <- dplyr::bind_rows(all_tracks)
  }

  distance_deg <- function(lon1, lat1, lon2, lat2) {
    sqrt((lon1 - lon2)^2 + (lat1 - lat2)^2)
  }

  track_list <- track_df |>
    dplyr::arrange(track_order, step) |>
    dplyr::group_split(track_order)

  ordered_tracks <- list()
  ordered_tracks[[1]] <- track_list[[1]]

  if (length(track_list) > 1) {

    for (i in 2:length(track_list)) {

      previous_track <- ordered_tracks[[i - 1]]
      current_track <- track_list[[i]]

      previous_end <- previous_track[nrow(previous_track), ]
      current_start <- current_track[1, ]
      current_end <- current_track[nrow(current_track), ]

      dist_to_start <- distance_deg(
        previous_end$lon,
        previous_end$lat,
        current_start$lon,
        current_start$lat
      )

      dist_to_end <- distance_deg(
        previous_end$lon,
        previous_end$lat,
        current_end$lon,
        current_end$lat
      )

      if (dist_to_end < dist_to_start) {
        current_track <- current_track |>
          dplyr::arrange(dplyr::desc(step)) |>
          dplyr::mutate(
            step = dplyr::row_number()
          )
      }

      ordered_tracks[[i]] <- current_track
    }
  }

  track_df <- dplyr::bind_rows(ordered_tracks) |>
    dplyr::mutate(
      global_step = dplyr::row_number()
    )

  track_json <- jsonlite::toJSON(
    track_df,
    dataframe = "rows",
    auto_unbox = TRUE
  )

  # ------------------------------------------------------------
  # HTML animation
  # ------------------------------------------------------------

  html_code <- paste0('
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>ICESat-2 Orbit Animation</title>

<style>
html, body {
  margin: 0;
  padding: 0;
  width: 100%;
  height: 100%;
  overflow: hidden;
  background: black;
  font-family: Arial, sans-serif;
}

canvas {
  display: block;
}

#controls {
  position: absolute;
  top: 15px;
  left: 15px;
  z-index: 10;
  background: rgba(0,0,0,0.75);
  color: white;
  padding: 14px;
  border-radius: 10px;
  width: 450px;
  font-size: 14px;
}

button {
  margin: 3px;
  padding: 6px 10px;
  border: none;
  border-radius: 5px;
  cursor: pointer;
}

input[type=range] {
  width: 330px;
}

#timeLabel,
#trackLabel,
#status {
  margin-top: 8px;
}
</style>
</head>

<body>

<div id="controls">

<b>ICESat-2 Orbit Animation</b>

<br><br>

<button onclick="playAnimation()">Play</button>
<button onclick="pauseAnimation()">Pause</button>
<button onclick="resetAnimation()">Reset</button>

<br><br>

Track speed:<br>
<input type="range" min="1" max="15" value="', track_speed, '" id="trackSpeed">
<span id="trackSpeedValue">', track_speed, '</span>

<br><br>

Earth rotation speed:<br>
<input type="range" min="0" max="20" value="', earth_rotation_speed, '" id="rotationSpeed">
<span id="rotationSpeedValue">', earth_rotation_speed, '</span>

<div id="trackLabel">Track:</div>
<div id="timeLabel">Time:</div>
<div id="status">Loading...</div>

</div>

<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>

<script>

const trackData = ', track_json, ';

const textureFileName = "', earth_texture, '";

// ------------------------------------------------------------
// Scene
// ------------------------------------------------------------

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x000000);

// ------------------------------------------------------------
// Camera
// ------------------------------------------------------------

const camera = new THREE.PerspectiveCamera(
  45,
  window.innerWidth / window.innerHeight,
  0.1,
  1000
);

camera.position.set(0, 0, 5);

// ------------------------------------------------------------
// Renderer
// ------------------------------------------------------------

const renderer = new THREE.WebGLRenderer({
  antialias: true
});

renderer.setSize(
  window.innerWidth,
  window.innerHeight
);

renderer.setPixelRatio(
  window.devicePixelRatio
);

document.body.appendChild(
  renderer.domElement
);

// ------------------------------------------------------------
// Controls
// ------------------------------------------------------------

const controls = new THREE.OrbitControls(
  camera,
  renderer.domElement
);

controls.enableDamping = true;
controls.dampingFactor = 0.05;
controls.enableZoom = true;

// ------------------------------------------------------------
// Lights
// ------------------------------------------------------------

const ambientLight = new THREE.AmbientLight(
  0xffffff,
  2.0
);

scene.add(ambientLight);

// ------------------------------------------------------------
// Earth group
// ------------------------------------------------------------

const earthGroup = new THREE.Group();
scene.add(earthGroup);

// ------------------------------------------------------------
// Earth globe
// ------------------------------------------------------------

const earthRadius = 1.5;
const orbitRadius = earthRadius + 0.35;

const earthGeometry = new THREE.SphereGeometry(
  earthRadius,
  256,
  256
);

const fallbackMaterial = new THREE.MeshBasicMaterial({
  color: 0x1e66b1
});

const earth = new THREE.Mesh(
  earthGeometry,
  fallbackMaterial
);

earthGroup.add(earth);

// ------------------------------------------------------------
// Earth texture
// ------------------------------------------------------------

new THREE.TextureLoader().load(
  "./" + textureFileName,
  function(texture) {

    texture.minFilter = THREE.LinearFilter;
    texture.magFilter = THREE.LinearFilter;
    texture.generateMipmaps = false;
    texture.needsUpdate = true;

    earth.material.dispose();

    earth.material = new THREE.MeshBasicMaterial({
      map: texture
    });

    document.getElementById("status").innerHTML =
      "Earth texture loaded";
  }
);

// ------------------------------------------------------------
// Stars
// ------------------------------------------------------------

const starGeometry = new THREE.BufferGeometry();
const starPositions = [];

for (let i = 0; i < 7000; i++) {

  starPositions.push(
    (Math.random() - 0.5) * 120,
    (Math.random() - 0.5) * 120,
    (Math.random() - 0.5) * 120
  );
}

starGeometry.setAttribute(
  "position",
  new THREE.Float32BufferAttribute(
    starPositions,
    3
  )
);

const starMaterial = new THREE.PointsMaterial({
  color: 0xffffff,
  size: 0.04
});

scene.add(
  new THREE.Points(
    starGeometry,
    starMaterial
  )
);

// ------------------------------------------------------------
// Coordinate conversion
// ------------------------------------------------------------

function latLonToVector3(lat, lon, radius) {

  const phi = (90 - lat) * Math.PI / 180;
  const theta = (lon + 180) * Math.PI / 180;

  const x = -radius * Math.sin(phi) * Math.cos(theta);
  const z =  radius * Math.sin(phi) * Math.sin(theta);
  const y =  radius * Math.cos(phi);

  return new THREE.Vector3(
    x,
    y,
    z
  );
}

// ------------------------------------------------------------
// ICESat-2 white sphere
// ------------------------------------------------------------

const satelliteGeometry = new THREE.SphereGeometry(
  0.08,
  48,
  48
);

const satelliteMaterial = new THREE.MeshBasicMaterial({
  color: 0xffffff
});

const satellite = new THREE.Mesh(
  satelliteGeometry,
  satelliteMaterial
);

earthGroup.add(satellite);

// ------------------------------------------------------------
// ICESat-2 text label above the sphere
// ------------------------------------------------------------

const labelCanvas = document.createElement("canvas");
labelCanvas.width = 1024;
labelCanvas.height = 256;

const labelContext = labelCanvas.getContext("2d");

labelContext.clearRect(
  0,
  0,
  labelCanvas.width,
  labelCanvas.height
);

labelContext.font = "bold 90px Arial";
labelContext.textAlign = "center";
labelContext.textBaseline = "middle";

labelContext.lineWidth = 10;
labelContext.strokeStyle = "black";

labelContext.strokeText(
  "ICESat-2",
  labelCanvas.width / 2,
  labelCanvas.height / 2
);

labelContext.fillStyle = "white";

labelContext.fillText(
  "ICESat-2",
  labelCanvas.width / 2,
  labelCanvas.height / 2
);

const labelTexture = new THREE.CanvasTexture(labelCanvas);

const labelMaterial = new THREE.SpriteMaterial({
  map: labelTexture,
  transparent: true,
  depthWrite: false,
  depthTest: false
});

const satelliteLabel = new THREE.Sprite(
  labelMaterial
);

satelliteLabel.scale.set(
  0.65,
  0.16,
  1
);

earthGroup.add(satelliteLabel);

// ------------------------------------------------------------
// Laser beam pointing from satellite to Earth
// ------------------------------------------------------------

let laserBeam = null;

function updateLaserBeam(satPos) {

  if (laserBeam) {

    earthGroup.remove(laserBeam);
    laserBeam.geometry.dispose();
    laserBeam.material.dispose();
    laserBeam = null;
  }

  const groundPoint = satPos.clone()
    .normalize()
    .multiplyScalar(
      earthRadius + 0.01
    );

  const geometry = new THREE.BufferGeometry()
    .setFromPoints([
      satPos,
      groundPoint
    ]);

  const material = new THREE.LineBasicMaterial({
    color: 0x00ff00,
    transparent: true,
    opacity: 0.95
  });

  laserBeam = new THREE.Line(
    geometry,
    material
  );

  earthGroup.add(laserBeam);
}

// ------------------------------------------------------------
// Track animation
// ------------------------------------------------------------

let pointIndex = 0;
let running = false;

let trackPoints = [];
let activeTrackLine = null;

const trackColor = 0x00ff00;

function clearTrackLine() {

  if (activeTrackLine) {

    earthGroup.remove(activeTrackLine);
    activeTrackLine.geometry.dispose();
    activeTrackLine.material.dispose();
    activeTrackLine = null;
  }
}

function updateTrackLine() {

  clearTrackLine();

  if (trackPoints.length < 2) return;

  const geometry = new THREE.BufferGeometry()
    .setFromPoints(trackPoints);

  const material = new THREE.LineBasicMaterial({
    color: trackColor
  });

  activeTrackLine = new THREE.Line(
    geometry,
    material
  );

  earthGroup.add(activeTrackLine);
}

function updateSatellitePosition(pos) {

  satellite.position.copy(pos);

  const labelOffset = pos.clone()
    .normalize()
    .multiplyScalar(0.22);

  satelliteLabel.position.copy(
    pos.clone().add(labelOffset)
  );

  updateLaserBeam(pos);
}

function addCurrentPoint() {

  if (pointIndex >= trackData.length) {

    running = false;

    document.getElementById("status").innerHTML =
      "Track completed";

    return;
  }

  const p = trackData[pointIndex];

  const pos = latLonToVector3(
    p.lat,
    p.lon,
    orbitRadius
  );

  trackPoints.push(pos);

  updateSatellitePosition(pos);

  document.getElementById("trackLabel").innerHTML =
    "Point: " +
    (pointIndex + 1) +
    " of " +
    trackData.length +
    " | Track: " +
    p.track_id;

  document.getElementById("timeLabel").innerHTML =
    "Time: " +
    p.time_label;

  pointIndex++;

  updateTrackLine();
}

function playAnimation() {

  running = true;

  document.getElementById("status").innerHTML =
    "Animation running";
}

function pauseAnimation() {

  running = false;

  document.getElementById("status").innerHTML =
    "Animation paused";
}

function resetAnimation() {

  running = false;
  pointIndex = 0;
  trackPoints = [];

  clearTrackLine();

  if (laserBeam) {

    earthGroup.remove(laserBeam);
    laserBeam.geometry.dispose();
    laserBeam.material.dispose();
    laserBeam = null;
  }

  if (trackData.length > 0) {
    addCurrentPoint();
  }

  document.getElementById("status").innerHTML =
    "Animation reset";
}

document.getElementById("trackSpeed")
  .addEventListener("input", function() {

    document.getElementById("trackSpeedValue").innerHTML =
      this.value;
});

document.getElementById("rotationSpeed")
  .addEventListener("input", function() {

    document.getElementById("rotationSpeedValue").innerHTML =
      this.value;
});

if (trackData.length > 0) {

  addCurrentPoint();

} else {

  document.getElementById("status").innerHTML =
    "No points found";
}

// ------------------------------------------------------------
// Animation loop
// ------------------------------------------------------------

function animate() {

  requestAnimationFrame(animate);

  const rotationSpeed = Number(
    document.getElementById("rotationSpeed").value
  );

  earthGroup.rotation.y +=
    0.0008 * rotationSpeed;

  if (
    running &&
    pointIndex < trackData.length
  ) {

    const trackSpeed = Number(
      document.getElementById("trackSpeed").value
    );

    for (
      let s = 0;
      s < trackSpeed;
      s++
    ) {

      if (
        running &&
        pointIndex < trackData.length
      ) {

        addCurrentPoint();
      }
    }
  }

  controls.update();

  renderer.render(
    scene,
    camera
  );
}

animate();

window.addEventListener(
  "resize",
  function() {

    camera.aspect =
      window.innerWidth /
      window.innerHeight;

    camera.updateProjectionMatrix();

    renderer.setSize(
      window.innerWidth,
      window.innerHeight
    );
});

</script>
</body>
</html>
')

  out_file <- file.path(
    output_dir,
    output_file
  )

  writeLines(
    html_code,
    out_file
  )

  if (launch) {
    try(
      servr::daemon_stop(),
      silent = TRUE
    )

    servr::httd(
      dir = output_dir,
      port = 6733
    )

    utils::browseURL(
      paste0(
        "http://127.0.0.1:6733/",
        output_file
      )
    )
  }

  invisible(out_file)
}

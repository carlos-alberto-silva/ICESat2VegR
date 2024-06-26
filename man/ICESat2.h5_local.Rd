% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.icesat2.h5_local.R
\name{ICESat2.h5_local}
\alias{ICESat2.h5_local}
\title{The class representing the h5 file opened from local files.}
\description{
The class representing the h5 file opened from local files.

The class representing the h5 file opened from local files.
}
\details{
The variants \verb{_cloud} and \verb{_local} allows all the other functions
to use generic calls using the same interface, with each class
implementation is provided accordingly.

The regular usage does not require the user to work with those classes
as most other provided functions will actually give access to the most
common necessities for working with ICESat-2 data.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{h5}}{A pointer to \code{\link[hdf5r:H5File-class]{hdf5r::H5File}} in case
the user wants to access features not implemented yet}

\item{\code{beams}}{The \code{\linkS4class{character}} vector of beams available for the granule.}

\item{\code{strong_beams}}{The \code{\linkS4class{character}} vector of strong beams calculated using orbit_info}

\item{\code{weak_beams}}{The \code{\linkS4class{character}} vector of weak beams calculated using orbit_info}

\item{\code{isOpen}}{A flag to indicate if the file pointer has already been closed}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-ICESat2.h5_local-new}{\code{ICESat2.h5_local$new()}}
\item \href{#method-ICESat2.h5_local-ls}{\code{ICESat2.h5_local$ls()}}
\item \href{#method-ICESat2.h5_local-ls_groups}{\code{ICESat2.h5_local$ls_groups()}}
\item \href{#method-ICESat2.h5_local-ls_attrs}{\code{ICESat2.h5_local$ls_attrs()}}
\item \href{#method-ICESat2.h5_local-dt_datasets}{\code{ICESat2.h5_local$dt_datasets()}}
\item \href{#method-ICESat2.h5_local-exists}{\code{ICESat2.h5_local$exists()}}
\item \href{#method-ICESat2.h5_local-attr}{\code{ICESat2.h5_local$attr()}}
\item \href{#method-ICESat2.h5_local-close_all}{\code{ICESat2.h5_local$close_all()}}
\item \href{#method-ICESat2.h5_local-print}{\code{ICESat2.h5_local$print()}}
\item \href{#method-ICESat2.h5_local-clone}{\code{ICESat2.h5_local$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-new"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-new}{}}}
\subsection{Method \code{new()}}{
Direct initialization should not be used, it is handled by \code{ATL0X_read()}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$new(h5)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{h5}}{the result of the \code{\link[=ATLAS_dataFinder]{ATLAS_dataFinder()}}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The class object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-ls"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-ls}{}}}
\subsection{Method \code{ls()}}{
Lists the groups and datasets that are within current group
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$ls()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
List the groups and datasets within the current path
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-ls_groups"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-ls_groups}{}}}
\subsection{Method \code{ls_groups()}}{
Lists all grouops recursively
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$ls_groups(recursive = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{recursive}}{\code{\linkS4class{logical}}, default FALSE. If TRUE it will list
groups recursively and return the full path.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The \code{\linkS4class{character}} representing
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-ls_attrs"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-ls_attrs}{}}}
\subsection{Method \code{ls_attrs()}}{
Lists the available attributes
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$ls_attrs()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{\linkS4class{character}} vector of attributes available
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-dt_datasets"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-dt_datasets}{}}}
\subsection{Method \code{dt_datasets()}}{
Get datasets as data.table with columns (name, dataset.dims, rank)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$dt_datasets(recursive = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{recursive}}{\code{\linkS4class{logical}}, default FALSE. If TRUE recursively searches
and returns the full path.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{\link[data.table:data.table]{data.table::data.table}} with the columns (name, dataset.dims, rank)
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-exists"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-exists}{}}}
\subsection{Method \code{exists()}}{
Checks if a supplied group/dataset exist within the H5
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$exists(path)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{\code{\linkS4class{character}} with the path to test}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{\linkS4class{logical}} TRUE or FALSE if it exists or not
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-attr"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-attr}{}}}
\subsection{Method \code{attr()}}{
Read an attribute from h5
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$attr(attribute)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{attribute}}{\code{\linkS4class{character}} the address of the attribute to open}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-close_all"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-close_all}{}}}
\subsection{Method \code{close_all()}}{
Safely closes the h5 file pointer
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$close_all(silent = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{silent}}{\code{\linkS4class{logical}}, default TRUE. Will cast warning messages if silent = FALSE.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Nothing, just closes the file
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-print"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-print}{}}}
\subsection{Method \code{print()}}{
Prints the data in a friendly manner
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Added for compatibility purposes with generic signature}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Outputs information about object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_local-clone"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_local-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_local$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

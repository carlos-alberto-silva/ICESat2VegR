% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.icesat2.h5_cloud.R
\name{ICESat2.h5_cloud}
\alias{ICESat2.h5_cloud}
\title{The class representing the h5 file opened from the cloud for cloud computing.}
\description{
The class representing the h5 file opened from the cloud for cloud computing.

The class representing the h5 file opened from the cloud for cloud computing.
}
\details{
Besides representing h5 files, it is also used to represent groups
opened with \verb{[[]]}.

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
\item{\code{h5}}{A pointer to \code{h5py} in case the user wants to access features not
implemented yet}

\item{\code{beams}}{The \code{\linkS4class{character}} vector of beams available for the granule.}

\item{\code{strong_beams}}{The \code{\linkS4class{character}} vector of strong beams calculated using orbit_info}

\item{\code{weak_beams}}{The \code{\linkS4class{character}} vector of weak beams calculated using orbit_info}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-ICESat2.h5_cloud-new}{\code{ICESat2.h5_cloud$new()}}
\item \href{#method-ICESat2.h5_cloud-ls}{\code{ICESat2.h5_cloud$ls()}}
\item \href{#method-ICESat2.h5_cloud-ls_groups}{\code{ICESat2.h5_cloud$ls_groups()}}
\item \href{#method-ICESat2.h5_cloud-ls_attrs}{\code{ICESat2.h5_cloud$ls_attrs()}}
\item \href{#method-ICESat2.h5_cloud-dt_datasets}{\code{ICESat2.h5_cloud$dt_datasets()}}
\item \href{#method-ICESat2.h5_cloud-exists}{\code{ICESat2.h5_cloud$exists()}}
\item \href{#method-ICESat2.h5_cloud-attr}{\code{ICESat2.h5_cloud$attr()}}
\item \href{#method-ICESat2.h5_cloud-close_all}{\code{ICESat2.h5_cloud$close_all()}}
\item \href{#method-ICESat2.h5_cloud-print}{\code{ICESat2.h5_cloud$print()}}
\item \href{#method-ICESat2.h5_cloud-clone}{\code{ICESat2.h5_cloud$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_cloud-new"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-new}{}}}
\subsection{Method \code{new()}}{
Direct initialization should not be used, it is handled by \code{ATL0X_read()}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$new(h5)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{h5}}{the result of the \code{\link[=ATLAS_dataFinder]{ATLAS_dataFinder()}} with \code{cloud_computing = TRUE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The class object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_cloud-ls"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-ls}{}}}
\subsection{Method \code{ls()}}{
Lists the groups and datasets that are within current group
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$ls()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
List the groups and datasets within the current path
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_cloud-ls_groups"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-ls_groups}{}}}
\subsection{Method \code{ls_groups()}}{
Lists all grouops recursively
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$ls_groups(recursive = FALSE)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-ICESat2.h5_cloud-ls_attrs"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-ls_attrs}{}}}
\subsection{Method \code{ls_attrs()}}{
Lists the available attributes
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$ls_attrs()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{\linkS4class{character}} vector of attributes available
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_cloud-dt_datasets"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-dt_datasets}{}}}
\subsection{Method \code{dt_datasets()}}{
Get datasets as data.table with columns (name, dataset.dims, rank)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$dt_datasets(recursive = FALSE)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-ICESat2.h5_cloud-exists"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-exists}{}}}
\subsection{Method \code{exists()}}{
Checks if a supplied group/dataset exist within the H5
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$exists(path)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-ICESat2.h5_cloud-attr"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-attr}{}}}
\subsection{Method \code{attr()}}{
Read an attribute from h5
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$attr(attribute)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-ICESat2.h5_cloud-close_all"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-close_all}{}}}
\subsection{Method \code{close_all()}}{
Safely closes the h5 file pointer
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$close_all()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Nothing, just closes the file
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ICESat2.h5_cloud-print"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-print}{}}}
\subsection{Method \code{print()}}{
Prints the data in a friendly manner
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$print(...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-ICESat2.h5_cloud-clone"></a>}}
\if{latex}{\out{\hypertarget{method-ICESat2.h5_cloud-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ICESat2.h5_cloud$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

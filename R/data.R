#' Davis Berlind's Daily Step Count
#'
#' Davis Berlind's total daily steps between 9-12-2016 and 08-15-2025
#'
#' @format ## `daily_steps`
#' A data frame with 3,247 rows and 2 columns:
#' \describe{
#'   \item{day}{Date}
#'   \item{total_steps}{Total steps in day}
#' }
#' @source My Apple Health app.
"daily_steps"

#' Ion Channel Voltage
#'
#' Ionic current recordings from a porin in the outer membrane of a bacterial
#' cell using the voltage clamp technique. Data have been processed so that they
#' include every 11th observation.
#'
#' @format ## `ion_channel`
#' A data frame with 2,956 rows and 1 columns:
#' \describe{
#'   \item{y}{Ionic current recordings}
#' }
#' @source \href{https://www.uni-goettingen.de/en/213067.html}{Steinem lab (Institute of Organic and Biomolecular Chemistry, University of Gottingen)}
"ion_channel"

#' Shankle Oil Well Lithology
#'
#' Six petrophysical measures collected by the Kansas Geological Survey by
#' drilling into the B1-C formations of the Shankle well in the Panoma oil and
#' gas field of southwest Kansas.
#'
#' @format ## `well_log`
#' A data frame with 346 rows and 8 columns:
#' \describe{
#'   \item{Facies}{1: "Nonmarine Sandstone, 2: "Nonmarine Coarse Siltstone, 3: Nonmarine Fine Siltstone, 4: Marine Siltstone/Shale, 5: Mudstone, 6: Wackestone, 7: Dolomite, 8: Packstone-Grainstone,8: Packstone-Grainstone"}
#'   \item{Depth}{Well depth in feet}
#'   \item{GR}{Gamma ray emissions}
#'   \item{Rt}{Deep resistivity}
#'   \item{DeltaPhi}{Difference of neutron and density porosity}
#'   \item{AvgPhi}{Average neutron and density porosity}
#'   \item{PE}{Photoelectric factor}
#'   \item{MnM}{Marine/non-marine indicator}
#' }
#' @source Bohling, G. and M. Dubois (2003). An integrated application of neural
#'   network and markov chain techniques to the prediction of lithofacies from
#'   well logs: Kansas geological survey open-file report 2003-50, 6 p. Group 6.
"well_log"

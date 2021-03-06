#' @title CPS-TUS data 2014-2015
#' @description Data from 2014-15 CPS Tobacco Use Supplement (TUS; U.S. Department of Commerce and U.S. Census Bureau 2016),
#' employed adult daily smokers (Non-Hispanic White males only). Missing data are removed from the dataset.
#' Due to the CRAN limitation of the size (5MB) of an R package, only half of the observations remained in this internal dataset for the purpose of illustration.
#' @format A data frame with 3922 rows and 167 variables:
#' \describe{
#'   \item{\code{PRTAGE}}{Age}
#'   \item{\code{PESEX}}{Gender, 0=Male, 1=Female}
#'   \item{\code{repwgt0}}{Sample main weights}
#'   \item{\code{repwgt1}}{BRR replicate weights}
#'   \item{\code{repwgt2}}{BRR replicate weights}
#'   \item{\code{repwgt3}}{BRR replicate weights}
#'   \item{\code{repwgt4}}{BRR replicate weights}
#'   \item{\code{repwgt5}}{BRR replicate weights}
#'   \item{\code{repwgt6}}{BRR replicate weights}
#'   \item{\code{repwgt7}}{BRR replicate weights}
#'   \item{\code{repwgt8}}{BRR replicate weights}
#'   \item{\code{repwgt9}}{BRR replicate weights}
#'   \item{\code{repwgt10}}{BRR replicate weights}
#'   \item{\code{repwgt11}}{BRR replicate weights}
#'   \item{\code{repwgt12}}{BRR replicate weights}
#'   \item{\code{repwgt13}}{BRR replicate weights}
#'   \item{\code{repwgt14}}{BRR replicate weights}
#'   \item{\code{repwgt15}}{BRR replicate weights}
#'   \item{\code{repwgt16}}{BRR replicate weights}
#'   \item{\code{repwgt17}}{BRR replicate weights}
#'   \item{\code{repwgt18}}{BRR replicate weights}
#'   \item{\code{repwgt19}}{BRR replicate weights}
#'   \item{\code{repwgt20}}{BRR replicate weights}
#'   \item{\code{repwgt21}}{BRR replicate weights}
#'   \item{\code{repwgt22}}{BRR replicate weights}
#'   \item{\code{repwgt23}}{BRR replicate weights}
#'   \item{\code{repwgt24}}{BRR replicate weights}
#'   \item{\code{repwgt25}}{BRR replicate weights}
#'   \item{\code{repwgt26}}{BRR replicate weights}
#'   \item{\code{repwgt27}}{BRR replicate weights}
#'   \item{\code{repwgt28}}{BRR replicate weights}
#'   \item{\code{repwgt29}}{BRR replicate weights}
#'   \item{\code{repwgt30}}{BRR replicate weights}
#'   \item{\code{repwgt31}}{BRR replicate weights}
#'   \item{\code{repwgt32}}{BRR replicate weights}
#'   \item{\code{repwgt33}}{BRR replicate weights}
#'   \item{\code{repwgt34}}{BRR replicate weights}
#'   \item{\code{repwgt35}}{BRR replicate weights}
#'   \item{\code{repwgt36}}{BRR replicate weights}
#'   \item{\code{repwgt37}}{BRR replicate weights}
#'   \item{\code{repwgt38}}{BRR replicate weights}
#'   \item{\code{repwgt39}}{BRR replicate weights}
#'   \item{\code{repwgt40}}{BRR replicate weights}
#'   \item{\code{repwgt41}}{BRR replicate weights}
#'   \item{\code{repwgt42}}{BRR replicate weights}
#'   \item{\code{repwgt43}}{BRR replicate weights}
#'   \item{\code{repwgt44}}{BRR replicate weights}
#'   \item{\code{repwgt45}}{BRR replicate weights}
#'   \item{\code{repwgt46}}{BRR replicate weights}
#'   \item{\code{repwgt47}}{BRR replicate weights}
#'   \item{\code{repwgt48}}{BRR replicate weights}
#'   \item{\code{repwgt49}}{BRR replicate weights}
#'   \item{\code{repwgt50}}{BRR replicate weights}
#'   \item{\code{repwgt51}}{BRR replicate weights}
#'   \item{\code{repwgt52}}{BRR replicate weights}
#'   \item{\code{repwgt53}}{BRR replicate weights}
#'   \item{\code{repwgt54}}{BRR replicate weights}
#'   \item{\code{repwgt55}}{BRR replicate weights}
#'   \item{\code{repwgt56}}{BRR replicate weights}
#'   \item{\code{repwgt57}}{BRR replicate weights}
#'   \item{\code{repwgt58}}{BRR replicate weights}
#'   \item{\code{repwgt59}}{BRR replicate weights}
#'   \item{\code{repwgt60}}{BRR replicate weights}
#'   \item{\code{repwgt61}}{BRR replicate weights}
#'   \item{\code{repwgt62}}{BRR replicate weights}
#'   \item{\code{repwgt63}}{BRR replicate weights}
#'   \item{\code{repwgt64}}{BRR replicate weights}
#'   \item{\code{repwgt65}}{BRR replicate weights}
#'   \item{\code{repwgt66}}{BRR replicate weights}
#'   \item{\code{repwgt67}}{BRR replicate weights}
#'   \item{\code{repwgt68}}{BRR replicate weights}
#'   \item{\code{repwgt69}}{BRR replicate weights}
#'   \item{\code{repwgt70}}{BRR replicate weights}
#'   \item{\code{repwgt71}}{BRR replicate weights}
#'   \item{\code{repwgt72}}{BRR replicate weights}
#'   \item{\code{repwgt73}}{BRR replicate weights}
#'   \item{\code{repwgt74}}{BRR replicate weights}
#'   \item{\code{repwgt75}}{BRR replicate weights}
#'   \item{\code{repwgt76}}{BRR replicate weights}
#'   \item{\code{repwgt77}}{BRR replicate weights}
#'   \item{\code{repwgt78}}{BRR replicate weights}
#'   \item{\code{repwgt79}}{BRR replicate weights}
#'   \item{\code{repwgt80}}{BRR replicate weights}
#'   \item{\code{repwgt81}}{BRR replicate weights}
#'   \item{\code{repwgt82}}{BRR replicate weights}
#'   \item{\code{repwgt83}}{BRR replicate weights}
#'   \item{\code{repwgt84}}{BRR replicate weights}
#'   \item{\code{repwgt85}}{BRR replicate weights}
#'   \item{\code{repwgt86}}{BRR replicate weights}
#'   \item{\code{repwgt87}}{BRR replicate weights}
#'   \item{\code{repwgt88}}{BRR replicate weights}
#'   \item{\code{repwgt89}}{BRR replicate weights}
#'   \item{\code{repwgt90}}{BRR replicate weights}
#'   \item{\code{repwgt91}}{BRR replicate weights}
#'   \item{\code{repwgt92}}{BRR replicate weights}
#'   \item{\code{repwgt93}}{BRR replicate weights}
#'   \item{\code{repwgt94}}{BRR replicate weights}
#'   \item{\code{repwgt95}}{BRR replicate weights}
#'   \item{\code{repwgt96}}{BRR replicate weights}
#'   \item{\code{repwgt97}}{BRR replicate weights}
#'   \item{\code{repwgt98}}{BRR replicate weights}
#'   \item{\code{repwgt99}}{BRR replicate weights}
#'   \item{\code{repwgt100}}{BRR replicate weights}
#'   \item{\code{repwgt101}}{BRR replicate weights}
#'   \item{\code{repwgt102}}{BRR replicate weights}
#'   \item{\code{repwgt103}}{BRR replicate weights}
#'   \item{\code{repwgt104}}{BRR replicate weights}
#'   \item{\code{repwgt105}}{BRR replicate weights}
#'   \item{\code{repwgt106}}{BRR replicate weights}
#'   \item{\code{repwgt107}}{BRR replicate weights}
#'   \item{\code{repwgt108}}{BRR replicate weights}
#'   \item{\code{repwgt109}}{BRR replicate weights}
#'   \item{\code{repwgt110}}{BRR replicate weights}
#'   \item{\code{repwgt111}}{BRR replicate weights}
#'   \item{\code{repwgt112}}{BRR replicate weights}
#'   \item{\code{repwgt113}}{BRR replicate weights}
#'   \item{\code{repwgt114}}{BRR replicate weights}
#'   \item{\code{repwgt115}}{BRR replicate weights}
#'   \item{\code{repwgt116}}{BRR replicate weights}
#'   \item{\code{repwgt117}}{BRR replicate weights}
#'   \item{\code{repwgt118}}{BRR replicate weights}
#'   \item{\code{repwgt119}}{BRR replicate weights}
#'   \item{\code{repwgt120}}{BRR replicate weights}
#'   \item{\code{repwgt121}}{BRR replicate weights}
#'   \item{\code{repwgt122}}{BRR replicate weights}
#'   \item{\code{repwgt123}}{BRR replicate weights}
#'   \item{\code{repwgt124}}{BRR replicate weights}
#'   \item{\code{repwgt125}}{BRR replicate weights}
#'   \item{\code{repwgt126}}{BRR replicate weights}
#'   \item{\code{repwgt127}}{BRR replicate weights}
#'   \item{\code{repwgt128}}{BRR replicate weights}
#'   \item{\code{repwgt129}}{BRR replicate weights}
#'   \item{\code{repwgt130}}{BRR replicate weights}
#'   \item{\code{repwgt131}}{BRR replicate weights}
#'   \item{\code{repwgt132}}{BRR replicate weights}
#'   \item{\code{repwgt133}}{BRR replicate weights}
#'   \item{\code{repwgt134}}{BRR replicate weights}
#'   \item{\code{repwgt135}}{BRR replicate weights}
#'   \item{\code{repwgt136}}{BRR replicate weights}
#'   \item{\code{repwgt137}}{BRR replicate weights}
#'   \item{\code{repwgt138}}{BRR replicate weights}
#'   \item{\code{repwgt139}}{BRR replicate weights}
#'   \item{\code{repwgt140}}{BRR replicate weights}
#'   \item{\code{repwgt141}}{BRR replicate weights}
#'   \item{\code{repwgt142}}{BRR replicate weights}
#'   \item{\code{repwgt143}}{BRR replicate weights}
#'   \item{\code{repwgt144}}{BRR replicate weights}
#'   \item{\code{repwgt145}}{BRR replicate weights}
#'   \item{\code{repwgt146}}{BRR replicate weights}
#'   \item{\code{repwgt147}}{BRR replicate weights}
#'   \item{\code{repwgt148}}{BRR replicate weights}
#'   \item{\code{repwgt149}}{BRR replicate weights}
#'   \item{\code{repwgt150}}{BRR replicate weights}
#'   \item{\code{repwgt151}}{BRR replicate weights}
#'   \item{\code{repwgt152}}{BRR replicate weights}
#'   \item{\code{repwgt153}}{BRR replicate weights}
#'   \item{\code{repwgt154}}{BRR replicate weights}
#'   \item{\code{repwgt155}}{BRR replicate weights}
#'   \item{\code{repwgt156}}{BRR replicate weights}
#'   \item{\code{repwgt157}}{BRR replicate weights}
#'   \item{\code{repwgt158}}{BRR replicate weights}
#'   \item{\code{repwgt159}}{BRR replicate weights}
#'   \item{\code{repwgt160}}{BRR replicate weights}
#'   \item{\code{sp_adltban}}{Supporting toward smoking ban in adults-exclusive areas (e.g., casinos, bars)}
#'   \item{\code{sp_kidsban}}{Supporting toward smoking ban in kids-related areas (e.g., children playground)}
#'   \item{\code{workban}}{Score of smoking restriction score at workplace }
#'   \item{\code{numcg}}{Number of cigarretes smoked per day}
#'}
#' @source \url{https://cancercontrol.cancer.gov/brp/tcrb/tus-cps/questionnaires-data}
#' @references U.S. Department of Commerce, & U.S. Census Bureau. (2016). National Cancer Institute and Food and Drug Administration co-sponsored Tobacco Use Supplement to the Current Population Survey. 2014-15.
"MedData"

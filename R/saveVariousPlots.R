#' Saves ggplot plots with reasonable defaults in png and svg (png has explicit white background (transparency off))
#'
#' Figures look great when using cowplot theme as default.
#'
#' @param filename the filename base to save the figures
#' @param plotobj ggplot image object
#' @param scale default 1.5
#' @param width default 25
#' @param height default 12
#' @param units default cm
#' @param bg_png default png background white (alternative "transparent" or a color)
#' @param filetype override to custom ggsave filetype (default "")
#'
#' @returns ggplot figures saved in png (quick pastes) and svg (vector graphics) or alternative with custom ggplot2 filetype
#'
#' @export
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' 
#' @import ggplot2
#'
#' @note
#'   Updates:
#'          2025-02-18 Alx added default white background for PNG only (saved trasparent earlier).
#'          2024-12-21 Alx removed default pdf
#'          2020-03-02 Alx split bucket into functions
#'          2020-01-01 Alx initial checkin
#' 
# Alex Saves ggplot plots with reasonable defaults in png and svg using cowplot
SaveVariousPlots <- function(filename_base, plotobj, scale = 1.5, width = 40, height = 20, units = "cm" , bg_png="white", filetype="")
    {
    	require(ggplot2)
    	if( filetype > "" ) {
    		# only save in specified format
    		if( filetype == "png" | filetype == ".png" | filetype == ".PNG" ) {
	    		ggsave(
		            paste0(filename_base, filetype),
		            plotobj,
		            scale = scale,
		            width = width,
		            height = height,
		            units = units,
		            bg = bg_png,
		            limitsize = FALSE )
    		} else {
	    		ggsave(
		            paste0(filename_base, filetype),
		            plotobj,
		            scale = scale,
		            width = width,
		            height = height,
		            units = units,
		            limitsize = FALSE )
    		}
    	} else {
    		# save in the two default formats
	        # ggsave(
	        #     paste0(filename_base, ".pdf"),
	        #     plotobj,
	        #     scale = scale,
	        #     width = width,
	        #     height = height,
	        #     units = units,
	        #     limitsize = FALSE
	        # )
	        ggsave(
	            paste0(filename_base, ".png"),
	            plotobj,
	            scale = scale,
	            width = width,
	            height = height,
	            units = units,
	            bg = bg_png,
	            limitsize = FALSE
	        )
	        ggsave(
	            paste0(filename_base, ".svg"),
	            plotobj,
	            scale = scale,
	            width = width,
	            height = height,
	            units = units,
	            limitsize = FALSE
	        )
	    }
    }


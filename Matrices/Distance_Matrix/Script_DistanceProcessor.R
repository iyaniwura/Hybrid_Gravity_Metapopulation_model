# This file computes the distances between the LHA based on the data Michael sent to me

# Load ggplot2
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(extrafont)
library(ggthemes)
library(lemon)
library(latex2exp)
library(purrr)


data_path = "/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/Manuscript_Files/DistanceMatrix/New_Distance_File/"
output_path = "/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/Manuscript_Files/DistanceMatrix/New_Distance_File/"

RawDist_Data <- read.csv(paste(data_path, "LHA_distances_for_metapop_model.csv", sep=""))

List_LHA = c("Abbotsford", "Agassiz/Harrison", "Burnaby", "Chilliwack", "Delta", "Hope",
						 "Langley", "Maple Ridge/Pitt Meadows", "Mission" , 
						 "New Westminster", "South Surrey/White Rock", "Surrey", "Tri-Cities" )

n_region = length(List_LHA)
Distance_Matrix = matrix(0, nrow = n_region, ncol = n_region)

for (i in 1:n_region){ 
						LHA_kth  <- RawDist_Data %>%
													filter(LHA_1 == List_LHA[i]) %>%
													filter(LHA_2 %in% List_LHA)
					
						LHA_k = LHA_kth[order(LHA_kth$LHA_2),]
						
						
						LHA_k_Dist <- round(LHA_k$DIST/1000)
						
						# arranging the distances in a matrix
						if (i == 1){
									# Distance_Matrix[1,i] = 0
									Distance_Matrix[(i+1):n_region, i] = LHA_k_Dist
						} else if (i >  1 & i < n_region) {
									Distance_Matrix[1:(i-1), i] = LHA_k_Dist[1:i-1]
									Distance_Matrix[(i+1):n_region, i] = LHA_k_Dist[i:(n_region-1)]
							
						} else {
									Distance_Matrix[1:(n_region-1), i] = LHA_k_Dist
						}
							
							
}




D_Mat = unname(as.matrix(Distance_Matrix))
D_Vec = as.vector(t(D_Mat)) 

##########################################################################################################################

LHA_vec <- c("Abbotsford", "Agassiz/Harrison", "Burnaby", "Chilliwack", "Delta", "Hope",
						 "Langley", "Mission", "Maple Ridge/Pitt Meadows", "New Westminster", "South Surrey/White Rock", "Surrey", "Tri-Cities")

NumLHA = 13
Output <-	data.frame(DeviceCount = round(D_Vec,3)) %>%
	tibble() %>%
	mutate(degree = rep(1:NumLHA, NumLHA), age=rep(1:NumLHA,each = NumLHA)) %>%
	ggplot() +
	geom_tile(aes(x=degree,y=age,fill=DeviceCount)) + geom_text(aes(x=degree,y=age,label = round(DeviceCount, 5))) +
	scale_fill_gradient(low = "white", high = "steelblue2") +
	scale_x_discrete(name ="Region of origin", limits = LHA_vec )  +
	scale_y_discrete( name ="Destination region", limits = LHA_vec )  +
	# ggtitle("01-May-2020 to 30-Oct-2020") +
	theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1),
				plot.title = element_text(hjust = 0.5),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), 
				legend.position= "none",
				panel.background = element_rect(colour = "black", size=2, fill=NA)) +
	guides(fill = guide_colourbar(title = "Distance (Km)" ,label = TRUE))



Output


png(filename=paste(output_path,"DistanceMatrix_NotScaled_New.png", sep=""), width = 700, height = 600)
plot(Output)
dev.off()

# saving the data in a csv file
Distance_Matrix_df = data.frame(Distance_Matrix)
filepath = paste(output_path, "Unscaled_DistanceMatrix.csv", sep = "")
write.csv(Distance_Matrix_df, filepath, row.names = FALSE)




########################### calculating theprobability matrix ####################
k = 1
f_ij = 1/(1 + D_Mat)^k
pi_ij = f_ij/colSums(f_ij)

Pi_vec = as.vector(t(pi_ij))
NumLHA = 13
Output <-	data.frame(DeviceCount = round(Pi_vec,3)) %>%
	tibble() %>%
	mutate(degree = rep(1:NumLHA, NumLHA), age=rep(1:NumLHA,each = NumLHA)) %>%
	ggplot() +
	geom_tile(aes(x=degree,y=age,fill=DeviceCount)) + geom_text(aes(x=degree,y=age,label = round(DeviceCount, 5))) +
	scale_fill_gradient(low = "white", high = "steelblue2") +
	scale_x_discrete(name ="Region of origin", limits = LHA_vec )  +
	scale_y_discrete( name ="Destination region",limits = LHA_vec )  +
	# ggtitle("01-May-2020 to 30-Oct-2020") +
	theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1),
				plot.title = element_text(hjust = 0.5),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				legend.position= "none",
				panel.background = element_rect(colour = "black", size=2, fill=NA)) +
	guides(fill = guide_colourbar(title = "" ,label = TRUE))



Output


png(filename=paste(output_path,"Scaled_DistanceMatrix_New.png", sep=""), width = 700, height = 600)
plot(Output)
dev.off()

# saving the data in a csv file
pi_ij_df = data.frame(pi_ij)
filepath = paste(output_path, "Scaled_DistanceMatrix.csv", sep = "")
write.csv(pi_ij_df, filepath, row.names = FALSE)








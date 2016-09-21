## define class "segment"
## the S4-class "segment"
## corresponds to a segment in a RT-PCR measurement
## members are "num_timepoints", the number of distinct 
## measurements in this segment, "dt", the time between
## disting measurements,"num_rep_meas", describing how many
## measurements were performed per time step
## and "meas_point", which should be either "END" or "START",
## defining when the measurements were taken
RTPCRSegment <- setClass(
  "RTPCRSegment",
    slots=c(
      num_timepoints="numeric",
      dt="numeric",
      num_rep_meas="numeric",
      meas_point="character"
    ),
    prototype=list(
      num_timepoints=240,
      dt=30,
      num_rep_meas=3,
      meas_point="END"
    )
)

## S4-class "RTPCRMeasurement" describes
## a complete RT-measurement 
## "num_wells" is the number of (measured) wells/samples
## "well_names" is a list of characters describing each well/sample
## "num_channels" is the number of measured channels, by default 3
## "channel_names" is a list of characters describing each channel,
## e.g. list("Cy3","Cy5","Cy35")
## "num_segs" specifies how many segments the whole measurement consists of
## "seg_list" is a list of all segments in this measurement
RTPCRMeasurement <- setClass(
  "RTPCRMeasurement",
    slots=c(
        num_wells="numeric",
        well_names="list",
        num_channels="numeric",
        channel_names="list",
        num_segs="numeric",
        seg_list = "list"
    ),
    prototype=list(
      num_wells=4,
      well_names=list("24hb_Polym_3nM","24hb_Polym_6nM"),
      num_channels=3,
      channel_names=list("Cy5","Cy3","Cy35"),
      seg_list=c(),
      num_segs=3
    )
)

## InitMeasurementProtocol is given a list of specifications
## of the current RTPCRMeasurement, that might look like the following:
## list_oomp <- lisT(240,30,3,"END",240,60,3,"END",240,120,3,"END"),
## "num_wells" is the number of wells in this measurement
## and "well_names", a list of character describing the names of the wells/samples
## a typical call of "InitMeasurementProtocol" would look like this:
##
## list_ommp <- list(240,30,3,"END",240,60,3,"END",240,120,3,"END")
## nwells <- 4
## well_names <- list("24hb_6nM","24hb_4nM","24hb_3nM","24hb_2nM")
## mp <- InitMeasurementProtocol(list_oomp,nwells,well_names)
## !! The object mp basically contains a complete description of the RT-PCR Experiment" 
## and will subsequently be passed to "LoadRTPCRData()"
InitMeasurementProtocol <- function(lommp,num_wells=2,well_names=list("24hb_6nM","24hb_3nM"))
{
  llist <- length(lommp)
  if(llist%%4!=0)
  {
    print("error: mm-protocol - list length is not a multiple of 4....")
  }
  num_segs <- llist/4;
  meas_prot <- new("RTPCRMeasurement")
  meas_prot@num_segs <- num_segs
  meas_prot@num_wells <- num_wells
  meas_prot@well_names <- well_names
  list_segs <- list()
  
  for(i in 1:num_segs)
  {
    curr_seg <- new("RTPCRSegment")
    curr_seg@num_timepoints<-lommp[[4*(i-1)+1]]
    curr_seg@dt<-lommp[[4*(i-1)+2]]
    curr_seg@num_rep_meas <- lommp[[4*(i-1)+3]]
    curr_seg@meas_point <- lommp[[4*(i-1)+4]]
    list_segs <- c(list_segs,curr_seg)
  }
  meas_prot@seg_list <- list_segs
  return(meas_prot)
}

## LoadRTPCRData(filname,mp) 
LoadRTPCRData <- function(filename,mp)
{
  require("gdata")
  df <- read.xls(filename)
  
  loc_num_segs <- mp@num_segs
  ##start_time=0
  curr_rtpcr_seg <- mp@seg_list[[1]]
  
  curr_dt <- curr_rtpcr_seg@dt
  curr_num_time_pts <- curr_rtpcr_seg@num_timepoints
  curr_num_rep_meas <- curr_rtpcr_seg@num_rep_meas
  time_seq<-rep(seq(from=curr_dt,to=(curr_num_time_pts*curr_dt),by=curr_dt),each=3)  

  for(i in 2:loc_num_segs)
  {
    curr_rtpcr_seg <- mp@seg_list[[i]]
    curr_dt <- curr_rtpcr_seg@dt
    
    start_time <- time_seq[length(time_seq)]+curr_dt
      
    curr_num_time_pts <- curr_rtpcr_seg@num_timepoints
    
    curr_num_rep_meas <- curr_rtpcr_seg@num_rep_meas
    
    curr_end_time <- start_time+(curr_num_time_pts-1)*curr_dt
  temp_time_seq <- rep(seq(from=start_time,to=curr_end_time,by=curr_dt),each=3)
    time_seq <- c(time_seq,temp_time_seq)
  }
  ## calculate the total number of data points....
  tot_num_dpps <- length(time_seq)
  num_wells <- mp@num_wells
  num_channels <-mp@num_channels
  num_segments <- mp@num_segs
  channel_names <-mp@channel_names
  
  data_array <- array(0,dim=c(num_wells,tot_num_dpps,num_channels,2))
  
  ## determine which columns correspond to which channel..:
  col_indices <- c(2,5,8)

  curr_seg_start <- 2
  seg_starts <- array(0,dim=c(num_segments))
  offset_curr_well_points <- 0         
  for(index_segment in 1:num_segments)
  {
    seg_starts[index_segment] <- curr_seg_start
    curr_seg <- mp@seg_list[[index_segment]]
    points_pwell_curr_seg <- curr_seg@num_timepoints*curr_seg@num_rep_meas
    curr_seg_start <- curr_seg_start + num_wells*( points_pwell_curr_seg+2)
    
    for(index_well in 1:num_wells)
    {
      ## get start pos
      curr_well_start_pos <- seg_starts[index_segment]+(index_well-1)*(points_pwell_curr_seg+2)
     ## print("current well start pos: ")
      ##print(curr_well_start_pos)
      ##print(df[curr_well_start_pos,])
      for(index_channel in 1:num_channels)
      {
        for(j in 1:points_pwell_curr_seg)
        {
        ##  data_array[index_well,(offset_curr_well_points+1):(offset_curr_well_points+points_pwell_curr_seg),index_channel,1]<-df[curr_well_start_pos:curr_well_start_pos+points_pwell_curr_seg-1,col_indices[index_channel]+1]
        ##  data_array[index_well,(offset_curr_well_points+1):(offset_curr_well_points+points_pwell_curr_seg),index_channel,2]<-df[curr_well_start_pos:curr_well_start_pos+points_pwell_curr_seg-1,col_indices[index_channel]+2]
          data_array[index_well,offset_curr_well_points+j,index_channel,1]<-df[curr_well_start_pos+j-1,col_indices[index_channel]+1]
          data_array[index_well,offset_curr_well_points+j,index_channel,2]<-df[curr_well_start_pos+j-1,col_indices[index_channel]+2]
        } 
      }
    }
    offset_curr_well_points <- offset_curr_well_points+points_pwell_curr_seg
  }
  print(seg_starts)
  for(i in 1:length(seg_starts))
  {
    print(df[seg_starts[i],])
  }
  
  print(col_indices)
  print(tot_num_dpps)  
  df_return <- data.frame(time_seq)
  well_name_list <- mp@well_names
  channel_name_list <-mp@channel_names
  for(index_well in 1:num_wells)
  {
    for(index_channel in 1:num_channels)
    {
      curr_name <- paste(well_name_list[[index_well]],channel_name_list[[index_channel]],sep="_")
      curr_name2 <- paste(well_name_list[[index_well]],channel_name_list[[index_channel]],sep="_TEMP_")
      print(curr_name)
      df_return[curr_name] <- data_array[index_well,,index_channel,1]
      df_return[curr_name2] <- data_array[index_well,,index_channel,2]
    }
  }
  return(df_return)
}

## PlotIndividualTraces from 
PlotIndividualTraces <- function(df_result,mp,list_of_channels)
{
  well_names <- mp@well_names
  channel_names <- mp@channel_names
  
  for(i in 1:mp@num_wells)
  {
    for(j in 1:length(list_of_channels))
    {
      curr_name <- paste(well_names[[i]],list_of_channels[[j]],sep="_")
      print(curr_name)
      #plot.new()
     
      plot(df_result$time_seq,df_result[,curr_name],xlab="time [s]",ylab="intensity [a.u.]")
      text(curr_name)
    }
  }
}

PlotAllChannels <- function(df_result,mp)
{
  well_names <- mp@well_names
  channel_names <- mp@channel_names
  ## loop through wells and channels and plot...
  for(i in 1:mp@num_wells)
  {
    for(j in 1:mp@num_channels)
    {
      curr_name <- paste(well_names[[i]],channel_names[[j]],sep="_")
      if(j==1)
      {
        plot(df_result$time_seq,df_result[,curr_name],ylim=c(0,10000))
      }
      else
      {
        points(df_result$time_seq,df_result[,curr_name])
      }
    }
  }
}

## read measurement protocol from ".xls"-file:
## list(240,30,3,"END",240,60,3,"END",240,120,3,"END")
Get_List_Ommp_From_File <- function(filename)
{
  require("gdata")
  df_protocol <- read.xls(filename)
  print(df_protocol)
  list_ommp <- list()
  
  num_rows <- nrow(df_protocol)
  for(i in 1:num_rows)
  {
    list_ommp <- append(list_ommp,df_protocol[i,5])
    list_ommp <- append(list_ommp,df_protocol[i,3])
    list_ommp <- append(list_ommp,df_protocol[i,6])
    list_ommp <- append(list_ommp,toString(df_protocol[i,7]))
  }
  return(list_ommp)
}

Plot_Segment <- function(time_seq,intens_data,list_ommp,index_segment)
{
  print(index_segment)
  offset<-1;
  if(index_segment>1)
  {
    for(i in 1:(index_segment-1))
    {
      offset <- offset + list_ommp[[(i-1)*4+1]]*list_ommp[[(i-1)*4+3]]
    }
    first_index <- offset  
  }
  else
  {
    first_index <- 1
  }
  num_points <- list_ommp[[(index_segment-1)*4+1]]*list_ommp[[(index_segment-1)*4+3]]
  
  last_index <- offset+num_points-1
  plot(time_seq[first_index:last_index],intens_data[first_index:last_index])
  print(offset)
}
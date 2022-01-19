rm(list=ls()) #reset working space
graphics.off() #closing current or all graphical windows


subst_data = read.table("marked_duplicates.error_rate_by_read_position.txt", header = TRUE, sep = "", dec = ".")



ones = subst_data[subst_data$read_number == 1 ,]

twos = subst_data[subst_data$read_number == 2 ,]


sums = rowSums(ones[,6:17])

sums2 = rowSums(twos[,6:17])


sum_column=colSums(subst_data[,6:17])

for(i in seq(1,length(sum_column),by=1)){
  print(sum_column[i])
}


#load necessary libraries
library(ggplot2)
library(reshape2)

#create data frame 
df = data.frame(position=ones$position,
                a_to_c=log(ones$a_to_c_error_rate),
                a_to_g=log(ones$a_to_g_error_rate),
                a_to_t=log(ones$a_to_t_error_rate),
                c_to_a=log(ones$c_to_a_error_rate),
                c_to_g=log(ones$c_to_g_error_rate),
                c_to_t=log(ones$c_to_t_error_rate),
                g_to_a=log(ones$g_to_a_error_rate),
                g_to_c=log(ones$g_to_c_error_rate),
                g_to_t=log(ones$g_to_t_error_rate),
                t_to_a=log(ones$t_to_a_error_rate),
                t_to_c=log(ones$t_to_c_error_rate),
                t_to_g=log(ones$t_to_g_error_rate),
                sum = log(sums))

#melt data frame into long format
df <- melt(df ,  id.vars = 'position', variable.name = 'substitution')

#create line plot for each column in data frame
ggplot(df, aes(position, value)) +
  geom_line(aes(colour = substitution)) + ggtitle("Read 1")



df = data.frame(position=twos$position,
                a_to_c=log(twos$a_to_c_error_rate),
                a_to_g=log(twos$a_to_g_error_rate),
                a_to_t=log(twos$a_to_t_error_rate),
                c_to_a=log(twos$c_to_a_error_rate),
                c_to_g=log(twos$c_to_g_error_rate),
                c_to_t=log(twos$c_to_t_error_rate),
                g_to_a=log(twos$g_to_a_error_rate),
                g_to_c=log(twos$g_to_c_error_rate),
                g_to_t=log(twos$g_to_t_error_rate),
                t_to_a=log(twos$t_to_a_error_rate),
                t_to_c=log(twos$t_to_c_error_rate),
                t_to_g=log(twos$t_to_g_error_rate),
                error = log(twos$error_rate))

#melt data frame into long format
df <- melt(df ,  id.vars = 'position', variable.name = 'substitution')

#create line plot for each column in data frame
ggplot(df, aes(position, value)) +
  geom_line(aes(colour = substitution)) + ggtitle("Read 2")



df = data.frame(position=twos$position,
                read1=log(sums),
                read2=log(sums2))

#melt data frame into long format
df <- melt(df ,  id.vars = 'position', variable.name = 'substitutions')

#create line plot for each column in data frame
ggplot(df, aes(position, value)) +
  geom_line(aes(colour = substitutions)) + ggtitle("Reads")

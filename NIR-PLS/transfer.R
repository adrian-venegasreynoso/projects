generate_validation <- function(msc_reference)
{

        raw_data <- read.csv("STD_SETS.csv", header=TRUE)
        raw_validation <- read.csv("VALSET.csv", header=TRUE)


        #Master
        inst1 <- raw_data %>% filter(Used.instrument=="Instr_1") %>% select(-Used.instrument, -Sample.name)
        #Slaves calibration
        inst2 <- raw_data %>% filter(Used.instrument=="Instr_2") %>% select(-Used.instrument, -Sample.name)
        inst3 <- raw_data %>% filter(Used.instrument=="Instr_3") %>% select(-Used.instrument, -Sample.name)
        inst4 <- raw_data %>% filter(Used.instrument=="Instr_4") %>% select(-Used.instrument, -Sample.name)
        inst5 <- raw_data %>% filter(Used.instrument=="Instr_5") %>% select(-Used.instrument, -Sample.name)
        inst6 <- raw_data %>% filter(Used.instrument=="Instr_6") %>% select(-Used.instrument, -Sample.name)
        inst7 <- raw_data %>% filter(Used.instrument=="Instr_7") %>% select(-Used.instrument, -Sample.name)
        #Slaves validation
        val1 <- raw_validation %>% filter(Used.instrument=="Instr_1") %>% select(-Used.instrument, -Sample.name)
        val2 <- raw_validation %>% filter(Used.instrument=="Instr_2") %>% select(-Used.instrument, -Sample.name)
        val3 <- raw_validation %>% filter(Used.instrument=="Instr_3") %>% select(-Used.instrument, -Sample.name)
        val4 <- raw_validation %>% filter(Used.instrument=="Instr_4") %>% select(-Used.instrument, -Sample.name)
        val5 <- raw_validation %>% filter(Used.instrument=="Instr_5") %>% select(-Used.instrument, -Sample.name)
        val6 <- raw_validation %>% filter(Used.instrument=="Instr_6") %>% select(-Used.instrument, -Sample.name)
        val7 <- raw_validation %>% filter(Used.instrument=="Instr_7") %>% select(-Used.instrument, -Sample.name)

        inst1 <- apply(inst1, 2, as.numeric)
        inst1 <- snv(inst1)
        inst1 <- my_savgol(inst1, 21,2,2)

        inst2 <- apply(inst2, 2, as.numeric)
        inst2 <- snv(inst2)
        inst2 <- my_savgol(inst2, 21,2,2)

        inst3 <- apply(inst3, 2, as.numeric)
        inst3 <- snv(inst3)
        inst3 <- my_savgol(inst3, 21,2,2)

        inst4 <- apply(inst4, 2, as.numeric)
        inst4 <- snv(inst4)
        inst4 <- my_savgol(inst4,21,2,2)

        inst5 <- apply(inst5, 2, as.numeric)
        inst5 <- snv(inst5)
        inst5 <- my_savgol(inst5, 21,2,2)

        inst6 <- apply(inst6, 2, as.numeric)
        inst6 <- snv(inst6)
        inst6 <- my_savgol(inst6, 21,2,2)

        inst7 <- apply(inst7, 2, as.numeric)
        inst7 <- snv(inst7)
        inst7 <- my_savgol(inst7, 21,2,2)


        val1 <- unname(val1)           
        val1 <- apply(val1, 2, as.numeric)
        val1 <- snv(val1)
        val1 <- my_savgol(val1, 21,2,2)

        val2 <- apply(val2, 2, as.numeric)
        val2 <- snv(val2)
        val2 <- my_savgol(val2, 21,2,2)

        val3 <- apply(val3, 2, as.numeric)
        val3 <- snv(val3)
        val3 <- my_savgol(val3, 21,2,2)

        val4 <- apply(val4, 2, as.numeric)
        val4 <- snv(val4)
        val4 <- my_savgol(val4, 21,2,2)

        val5 <- apply(val5, 2, as.numeric)
        val5 <- snv(val5)
        val5 <- my_savgol(val5, 21,2,2)

        val6 <- apply(val6, 2, as.numeric)
        val6 <- snv(val6)
        val6 <- my_savgol(val6, 21,2,2)

        val7 <- apply(val7, 2, as.numeric)
        val7 <- snv(val7)
        val7 <- my_savgol(val7, 21,2,2)


        val2 <- transfer_cal(inst1, inst2, wavelength, val2)
        val3 <- transfer_cal(inst1, inst3, wavelength, val3)
        val4 <- transfer_cal(inst1, inst4, wavelength, val4)
        val5 <- transfer_cal(inst1, inst5, wavelength, val5)
        val6 <- transfer_cal(inst1, inst6, wavelength, val6)
        val7 <- transfer_cal(inst1, inst7, wavelength, val7)

        X <-rbind(val1,val2, val3, val4, val5, val6, val7)
#        X <- my_savgol(X, 21,2,2)
        #spectra_plot(X, wavelength)
        #write.csv(X, "validation_set.csv")
        return(X)
}

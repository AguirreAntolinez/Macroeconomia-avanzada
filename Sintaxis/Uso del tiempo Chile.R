library(tidyverse)
uso_tiempo_chile<-read.csv2("https://raw.githubusercontent.com/AguirreAntolinez/Macroeconomia-avanzada/main/Datos%20de%20Chile/BASE_USUARIO%20ENUT%202015.csv")

uso_tiempo_chile<-uso_tiempo_chile %>% mutate(
  o61_1_2= ifelse(is.na(o61_1_2),0,o61_1_2),
  o61_2_2= ifelse(is.na(o61_2_2),0,o61_2_2),
  o62_1_2= ifelse(is.na(o62_1_2),0,o62_1_2),
  o62_2_2= ifelse(is.na(o62_2_2),0,o62_2_2),
  s21_1_2= ifelse(is.na(s21_1_2),0,s21_1_2),
  s21_2_2= ifelse(is.na(s21_2_2),0,s21_2_2),
  s22_1_2= ifelse(is.na(s22_1_2),0,s22_1_2),
  s22_2_2= ifelse(is.na(s22_2_2),0,s22_2_2),
  s23_1_2= ifelse(is.na(s23_1_2),0,s23_1_2),
  s23_2_2= ifelse(is.na(s23_2_2),0,s23_2_2),
  m11_1_2= ifelse(is.na(m11_1_2),0,m11_1_2),
  m11_2_2= ifelse(is.na(m11_2_2),0,m11_2_2),
  m12_1_2= ifelse(is.na(m12_1_2),0,m12_1_2),
  m12_2_2= ifelse(is.na(m12_2_2),0,m12_2_2),
  
  tiempo_consumo= 
    o61_1_2 + o61_2_2 + 
    o62_1_2 + o62_2_2 + 
    s21_1_2 + s21_2_2 + 
    s22_1_2 + s22_2_2 + 
    s23_1_2 + s23_2_2,
  
  tiempo_trabajo=
    m11_1_2 + m11_2_2 + 
    m12_1_2 + m12_2_2,
  
  tiempo_consumo= ifelse(tiempo_consumo==0,NA,tiempo_consumo),
  tiempo_trabajo= ifelse(tiempo_trabajo==0,NA,tiempo_trabajo)
)

summary(uso_tiempo_chile$tiempo_trabajo)
summary(uso_tiempo_chile$wgt2)
#Para el trabajo filtramos los que tengan peso y que esten entre 18 años y 65 años
uso_tiempo_chile_2<-uso_tiempo_chile%>% filter(!is.na(wgt2),c14_1_1>17, c14_1_1<66) 
summary(uso_tiempo_chile_2$tiempo_trabajo)
#Tiempo dedicado al trabajo a la semana
weighted.mean(x = uso_tiempo_chile_2$tiempo_trabajo,
              w = uso_tiempo_chile_2$wgt2,na.rm = TRUE)

#Tiempo dedicado al consumo a la semana
uso_tiempo_chile_3<-uso_tiempo_chile%>% filter(!is.na(wgt2)) 
weighted.mean(x = uso_tiempo_chile_3$tiempo_consumo,
              w = uso_tiempo_chile_3$wgt2,na.rm = TRUE)

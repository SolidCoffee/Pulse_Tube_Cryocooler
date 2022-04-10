#include "NTC_driver.h"

int RTC_C=A1;
int RTC_H=A2;

void setup() {
  // put your setup code here, to run once:
 
Serial.begin(9600);
pinMode(INPUT,RTC_C);
pinMode(INPUT,RTC_H);
}

void loop() {
  // put your main code here, to run repeatedly:
RTCval=analogRead(RTC_C);
RTCval_H=analogRead(RTC_H);

RTDTemp(RTCval);
RTDTemp(RTCval_H);

delay(50);
}

void RTDTemp(int valRTC)
{
  while(1)
  {
    if(valRTC == 0)
    {
      Serial.println("TOO HOT or error");
      break;
    }
    else if(valRTC == 1023)
    {
      Serial.println("Error or too cold");
      break;
    }
  
    if(valRTC <= R_RTD[i])
    {
      int x=(valRTC-R_RTD[i])/(R_RTD[i-1]-R_RTD[i]);
      int y=T_RTD[i-1]-T_RTD[i];
     Temp_C=(x*y)+T_RTD[i];
     Temp_K=Temp_C+273;
     Temp_F=(Temp_C*(9.0/5.0))+32.0;
     //Serial.print("Temp C =");
     //Serial.println(Temp_C);
     //Serial.print("Temp F =");
     //Serial.println(Temp_F);
     Serial.println(Temp_K);
     i=0;
     break;
    }
    i++;
  }
}

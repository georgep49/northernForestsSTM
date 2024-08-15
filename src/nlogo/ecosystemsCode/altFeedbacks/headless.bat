REM headless batch 

REM java  -Xmx1024M -cp C:\PROGRA~2\NETLOG~1.1\NetLogo.jar org.nlogo.headless.Main --model J:\private\SUTXCA~L\MGT6RW~J\FRWAXG~R\A5JIZE~I\altFeedbacks.nlogo --experiment predation --threads 4 --table "predationExpt.csv"

REM java  -Xmx1024M -cp C:\PROGRA~2\NETLOG~1.1\NetLogo.jar org.nlogo.headless.Main --model J:\private\SUTXCA~L\MGT6RW~J\FRWAXG~R\A5JIZE~I\altFeedbacks.nlogo --experiment consumed --threads 4  --table "consumeExpt.csv"

REM java  -Xmx1024M -cp C:\PROGRA~2\NETLOG~1.1\NetLogo.jar org.nlogo.headless.Main --model J:\private\SUTXCA~L\MGT6RW~J\FRWAXG~R\A5JIZE~I\altFeedbacks.nlogo --experiment fire-noinvasion --threads 4 --table "ffNoInv.csv"

REM java  -Xmx1024M -cp C:\PROGRA~2\NETLOG~1.1\NetLogo.jar org.nlogo.headless.Main --model J:\private\SUTXCA~L\MGT6RW~J\FRWAXG~R\A5JIZE~I\altFeedbacks.nlogo --experiment fire-invasion --threads 4 --table "ffInv.csv"

REM java  -Xmx1024M -cp C:\PROGRA~2\NETLOG~1.1\NetLogo.jar org.nlogo.headless.Main --model J:\private\SUTXCA~L\MGT6RW~J\FRWAXG~R\A5JIZE~I\altFeedbacks.nlogo --experiment ldd-by-pred --threads 4 --table "predationLDDExpt.csv"

java  -Xmx1024M -cp C:\PROGRA~2\NETLOG~1.1\NetLogo.jar org.nlogo.headless.Main --model C:\Users\gper020\Desktop\SC\ALTFEE~1\altFeedbacks.nlogo --experiment fire-invasion-regen --threads 2 --table "ffInvRegen.csv"
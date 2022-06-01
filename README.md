# 全反射透鏡 (TIR Lens)

## 前言
全反射透鏡 (Total Internal Reflection Lens, TIR Lens) 的功能為 : 將發散的光收攏，使其光源的光線準直。因此，TIR Lens 又可以稱作為準直器 (Collimator)。其目的就是提高照明時的光效率。如果針對特定範圍進行打光時，TIR Lens 可以減少光源的發散角度，將能量集中在目標範圍內部，因而減少能量損耗，來達到使用少量的光源具有一定程度的光照度。

---
## 介紹
TIR Lens 設計可以分成兩個部分 : **中央準直表面 (Center Surface)**與**全反射表面 (TIR Surface)**。由於 TIR Lens 是二次光學元件，一般是安裝在發光二極體 (Light Emission Diode, LED) 上，使 LED 的光線準直。而 LED 本身就是一個具有高指向性的發光元件，因此在小光束角，TIR Lens 就以中央準直表面進行角度的修正，使原本發散的角度修正成平行光出射，剩下的大光束角就藉由全反射表面修正成平行光出射，如下圖 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/How_TIR_Lens_works.png)

在此先將之後說明時的名詞進行定義，下圖表示 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Intordution_of_TIR_lens.png)

圖片中的名詞可以對到 maincode.m 的註解以及參數名詞，使用者可以參照設計需求來調整參數。程式參數控制一欄由下圖表示 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Parameter_setting.png)

Center Surface 與 TIR Surface 的計算方法是相同，為了方便了解，我做了一個TIR Surface 的流程圖，如下圖 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/flowchart.png)

按照上述流程就可以得到 Center Surface 與 TIR Surface，此程式會把取樣的每一點的切線向量 (Tangent Vector) 與法向量 (Normal Vector) 標記出來，如下圖 : 

<img src="https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Ray_Tracing_of_TIR_Surface.png" width="900" hight="900">

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Ray_Tracing_of_Center_Surface.png)

同時，會將所有的取樣的點進行高階曲線的擬合 (Curve Fitting)，並且將擬合的曲線與取樣的點進行比較，計算兩者的差異值 (Error Value)，如下圖表示 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/TIR_Surface_Fitting_Curve_Error_Value.png)
![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/TIR_Surface_Fitting_Curve.png)

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Center_Surface_Fitting_Curve_Error_Value.png)
![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Center_Surface_Fitting_Curve.png)

在參數控制一欄有 torlerant = 0.001，代表擬合的曲線與取樣的的最大的差異不可大過 0.001 Unit。而 torlerant 是可以根據使用者調整，數值越小，擬合的曲線階數就會越高。程式會在 Matlab 的 Command Windows 顯示不同階數的 Error Value，並且會列出曲線函數的係數與常數，如下圖所示 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/command_windows.png)

之後會在對 Fitting 的曲線進行驗證，計算經過兩曲線 (Center Surface 與 TIR Surface) 的光線出射角度是否有垂直出射 TIR Lens，如下圖表示 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Incident_Angle_vs_Reflection_Angle_on_TIR_Surface.png)
![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Incident_Angle_vs_Reflection_Angle_on_Center_Surface.png)

可以看到出射角度都非常接近 90 度 (z 軸方向)。(若想要更完美可以提高取樣點，以及減少 torlerant 值。)

最後，在參數控制一欄有 Propagate_distance，是可以顯示出出射 TIR Lens 後，傳播 Propagate_distance 的距離之光追跡效果，並且會顯示該截面的分布，如下表示 : 

![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/TIR_Lens_Ray_Tracing.png)
![image](https://github.com/YiChenLai/TIR_Lens/blob/master/pic/Surface_Cut_at_Porpagated_10_Unit_Distance.png)

---
## 結語
此程式目前可以根據使用者的設計條件進行**點光源**的 TIR Lens 生成，也就是說，此 TIR Lens 是使用理想光源設計出來的，而實際上 LED 的發光晶片都是面光源，代表此程式生成出來的透鏡沒有辦法達到最好的效果。但是使用者可以根據此程式再自行搭配演算法來修正表面來達到更好的效果。

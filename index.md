<img src="https://img.shields.io/badge/language-c++-brightgreen" alt="Language" style="zoom:150%;" /> <img src="https://img.shields.io/github/license/changqj/Mean-Value-Coordinates-in-3D" alt="github license" style="zoom:150%;" />



<div align="center">
    <img width="500px" src="../main/fig/MVC3D.png">
   <h1>
       <a href="https://github.com/changqj/Mean-Value-Coordinates-in-3D">Mean Value Coordinates in 3D [c++]</a>
    </h1>
   <h3>  | <a href="https://www.sciencedirect.com/science/article/pii/S0167839605000725"> Paper link on Sciencedirect </a>  | <a href="https://www.mn.uio.no/math/english/people/aca/michaelf/papers/mv3d.pdf"> Pdf version link </a> | </h3>
</div>



> This project implemented the **Mean Value Coordinates in 3D** algorithm in `c++`



### Usage

Assuming you have cloned this project in your home directory `d:\mvc3d`, we will next compile this project with `cmake`.

```powershell
$ cd d:\mvc3d
$ mkdir bin
$ cd bin
$ cmake ..
```

Once you have completed these steps, you can find the `.sln` file in the `\bin` folder, which is the project file that is successfully compiled for you. 

By compiling the project file, you can find the generated `mvc3d.exe` file in the `debug` or `release` folder. To run `mvc3d.exe`, the command line syntax is

```sh
$ mvc3d model_path x y z output_path
```

- model_path:  the current version only supports reading `.off` triangular mesh models. *Other formats will be supported as soon*.
- x & y & z: [x,y,z] is the point to calculate the barycentric coordinates.
- output_path: the calculated barycentric coordinates will be output to this file in text format.



**If you encounter any problems in the process of using this project, please do not hesitate to contact [me](mailto:qingjun.cn@gmail.com) or submit your problem to [issues](https://github.com/changqj/Mean-Value-Coordinates-in-3D/issues). Thanks for using.**

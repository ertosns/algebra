# algebra: C++ library for linear algebra

C++ linear algebra library required for [ kinetics ] (https://github.com/ertosns/kinetics.git) based of the API defined by  [_Modern Robotics: Mechanics, Planning, and Control_](https://modernrobotics.org) (Kevin Lynch and Frank Park, Cambridge University Press 2017).


# to deploy the library run:

```console
user@name:~$ . ./install.sh
```

# to verify is working, you can try an example examples/exercise.cpp

```console
user@name:~$ g++  exercise.cpp -lalgebra -I/usr/local/lib -I /usr/include/eigen3 -lpthread
```

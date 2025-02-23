import numpy as np
from matplotlib import pyplot as plt


def rotate_vec(vector, n):
    rotation_matrices = {
        0: np.array([[1, 0], [0, 1]]),    # 0 degrees
        1: np.array([[0, -1], [1, 0]]),   # 90 degrees
        2: np.array([[-1, 0], [0, -1]]),  # 180 degrees
        3: np.array([[0, 1], [-1, 0]])    # 270 degrees
    }
    
    roated_vector = rotation_matrices[n] @ vector
    return roated_vector


def segment_point(p1,p2):
    v1 = np.array(np.array(p2)-np.array(p1))
    relative_length = np.linalg.norm(v1)/4
    v1 = v1/np.linalg.norm(v1)
    v2 = np.array([0.0,1.0])
    angle = np.arctan2(v1[1], v1[0])
    angle_degrees = np.degrees(angle)
    angle_degrees = angle_degrees % 360
    n = round(angle_degrees / 90) % 4
    vectors_list = [[relative_length, 0.0], 
                    [0.0, relative_length], 
                    [relative_length, 0.0], 
                    [0.0, -relative_length], 
                    [0.0, -relative_length], 
                    [relative_length, 0.0], 
                    [0.0, relative_length]]
    point = np.array(p1)
    points_list = [p1]
    for vectors in vectors_list:
        roatated_vec = rotate_vec(vectors,n)
        point = point + roatated_vec
        points_list.append(np.round(point).tolist())
    return points_list

def new_line(line):
    new_line = []
    for i in range(0,len(line)-1):
        segment_points = segment_point(line[i], line[i+1]) 
        new_line = new_line + segment_points
    new_line.append(line[-1])
    return new_line

def fract(l, line):
    for _i in range(0,l):
        line = new_line(line)
    return line

def plotting_line(l, size):
    line = [[0.0,0.0], [size,0.0]]
    for i in range(0,l+1):
        line = fract(i,line)
        x_cords = [i[0] for i in line]
        y_cords = [i[1] for i in line]
        plt.plot(x_cords,y_cords)
    plt.show()

def plotting_square(l,size):
    lines = [
        [[0.0, 0.0], [size, 0.0]],  
        [[size, 0.0], [size, size]],
        [[size, size], [0.0, size]],
        [[0.0, size], [0.0, 0.0]]   
    ]
    for i in range(0,l+1):
        total_points = []
        for line in lines:
            line = fract(i,line)
            total_points = total_points + line
        x_cords = [i[0] for i in total_points]
        y_cords = [i[1] for i in total_points]
        plt.plot(x_cords,y_cords)
    plt.show()


def generate_squre(l,size):
    lines = [
        [[0.0, 0.0], [size, 0.0]],  
        [[size, 0.0], [size, size]],
        [[size, size], [0.0, size]],
        [[0.0, size], [0.0, 0.0]]
        ]
    total_points = []
    for line in lines:
        line = fract(l,line)
        total_points = total_points + line
    return total_points

    


problem_points = [
    (20, 134), # left wall
    (268, 151), # back wall
    (384, 240), # right wall
]

ok_points = [
    (114, 26), # ceiling
    (313, 349), # short box - right side
    (246, 315), # short box - front side 
    (228, 270), # short box - top side
    (82, 388), # floor
    (146, 210), # tall box - front side
    (94, 175), # tall box - left side
]

commands = [
    f"mitsuba cornell-box/scene-training.xml -o whatever.exr -D x={x} -D y={y} -D pdfCount=-1"
    for x, y in (problem_points + ok_points)
]

print(" && ".join(commands))

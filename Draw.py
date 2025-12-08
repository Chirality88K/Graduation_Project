import os
import matplotlib.pyplot as plt

def plot_txt_file(txt_path, img_path, is_cur=True):
    x, y = [], []
    with open(txt_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                print("value error!!!")
                continue
            try:
                x_val = float(parts[0])
                y_val = 0.0
                color = "black"
                title = ""
                label = ""
                if(is_cur):
                    y_val = float(parts[1])
                    color = "blue"
                    title = "Curvature"
                    label = "curvature"
                else:
                    y_val = float(parts[2])
                    color = "red"
                    title = "Torsion"
                    label = "torsion"
            except ValueError:
                continue
            x.append(x_val)
            y.append(y_val)
    if x and y:
        plt.figure()
        plt.plot(x, y, marker='o',color=color)
        plt.xlabel('t')
        plt.ylabel(label)
        plt.title(title)
        plt.grid(True)
        plt.savefig(img_path)
        plt.close()

def recursive_traverse(root_dir):
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # 跳过有 .visited 文件的文件夹
        if '.visited' in filenames and dirpath != root_dir:
            dirnames.clear()  # 不再遍历该文件夹的子文件夹
            continue

        for filename in filenames:
            if filename.lower().endswith('.txt'):
                txt_path = os.path.join(dirpath, filename)
                img_name = os.path.splitext(filename)[0] + 'cur.png'
                img_path = os.path.join(dirpath, img_name)
                plot_txt_file(txt_path, img_path, True)
                img_name = os.path.splitext(filename)[0] + 'tor.png'
                img_path = os.path.join(dirpath, img_name)
                plot_txt_file(txt_path, img_path, False)

        # 创建 .visited 文件，根目录不创建
        if dirpath != root_dir and '.visited' not in filenames:
            visited_path = os.path.join(dirpath, '.visited')
            with open(visited_path, 'w') as f:
                f.write('')

if __name__ == '__main__':
    root_folder = 'build/Table/'  # 改成你需要递归遍历的文件夹路径
    recursive_traverse(root_folder)
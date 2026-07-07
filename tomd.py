import os
import re

def optimize_matlab_to_md(input_dir, output_file):
    md_content = ["# MATLAB Model Source Code\n\n"]
    
    # 匹配 MATLAB 常见的注释块（如连续的 % 或许可证信息）
    # 这里主要针对文件开头的长注释进行压缩
    header_comment_pattern = re.compile(r'^(\s*%+.*\n)+')
    
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".m"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                    
                    # 1. 移除文件开头的冗余注释块（通常是 License 或作者信息）
                    content = header_comment_pattern.sub('', content, count=1)
                    
                    # 2. 压缩多余的空行（连续两个以上的换行符改为一个）
                    content = re.sub(r'\n\s*\n', '\n', content)
                    
                    # 3. 移除每行末尾的空格
                    lines = [line.rstrip() for line in content.split('\n')]
                    content = '\n'.join(lines)
                    
                    # 添加到 Markdown 列表
                    md_content.append(f"## File: {file}\n")
                    md_content.append("```matlab\n")
                    md_content.append(content)
                    md_content.append("\n```\n\n")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.writelines(md_content)
    
    print(f"成功！已将所有代码合并至: {output_file}")

# 使用示例
input_folder = r'D:\Research\GitHub\1-D-River-Sediment-Model\Coupled'  # 替换为你的代码路径
output_md = os.path.join(input_folder, 'Current_PDE_model.md')
optimize_matlab_to_md(input_folder, output_md)
# print(f"✅ 完成！文件已保存至：{output_md}")
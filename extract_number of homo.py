from collections import Counter
import re


def read_fasta(filename):
    """读取fasta文件，返回拼接好的序列字符串（忽略描述行）"""
    seq = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            seq.append(line)
    return ''.join(seq).upper()  # 转换为大写确保一致性


def analyze_homopolymers(seq, min_length=4):
    """分析同聚物，返回计数和位置信息"""
    pattern = r'(A{4,}|T{4,}|C{4,}|G{4,})'
    matches = re.finditer(pattern, seq)  # 使用finditer获取匹配位置

    counts = Counter()
    base_totals = Counter()
    all_homopolymers = []  # 存储所有长度≥min_length的同聚物信息

    for match in matches:
        base = match.group()[0]
        length = len(match.group())
        start_pos = match.start() + 1  # 转换为1-based坐标
        end_pos = match.end()

        counts[(base, length)] += 1
        base_totals[base] += length

        # 记录所有同聚物信息
        all_homopolymers.append({
            'base': f"{length}{base}",
            'start': start_pos,
            'end': end_pos,
            'length': length,
            'sequence': match.group()
        })

    return counts, base_totals, all_homopolymers


if __name__ == '__main__':
    fasta_file = "/Users/yzhou799/Desktop/sequence-36.fasta"
    output_file = "/Users/yzhou799/Desktop/fungi.txt"

    seq = read_fasta(fasta_file)
    counts, base_totals, all_homos = analyze_homopolymers(seq)

    # 生成TXT文件
    with open(output_file, 'w') as f:
        # 写入表头
        f.write("base\tstart\tend\n")

        # 写入所有同聚物信息
        for homo in all_homos:
            f.write(f"{homo['base']}\t{homo['start']}\t{homo['end']}\n")

    print(f"同聚物信息已保存到: {output_file}")

    # 输出统计信息到控制台
    print("\n--- 同聚物统计（长度≥4） ---")
    for base in ['A', 'T', 'C', 'G']:
        for length in sorted({l for (b, l) in counts if b == base}):
            print(f"{length}{base} Count: {counts[(base, length)]}")

    print("\n--- 总碱基个数统计（连续4个及以上） ---")
    for base in ['A', 'T', 'C', 'G']:
        print(f"{base} total count: {base_totals[base]}")

    # 输出长度≥10的同聚物信息
    long_homos = [h for h in all_homos if h['length'] >= 10]
    if long_homos:
        print(f"\n--- 长同聚物位置（长度≥10，共{len(long_homos)}处） ---")
        for i, homo in enumerate(long_homos, 1):
            print(f"{i}. {homo['base']} (位置: {homo['start']}-{homo['end']})")
            print(f"   序列: {homo['sequence']}")
    else:
        print("\n未发现长度≥10的同聚物")
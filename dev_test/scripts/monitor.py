"""Run a command and sample thread count / CPU cores used of its process tree.

Usage: python3 scripts/monitor.py <log.tsv> -- <command ...>
Prints a summary (peak threads, peak and mean cores busy) when the command exits.
"""
import os
import subprocess
import sys
import time

INTERVAL = 1.0
TICKS = os.sysconf("SC_CLK_TCK")


def descendants(root):
    children = {}
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            with open(f"/proc/{pid}/stat") as f:
                ppid = int(f.read().rsplit(")", 1)[1].split()[1])
        except (OSError, IndexError, ValueError):
            continue
        children.setdefault(ppid, []).append(int(pid))
    out, stack = [], [root]
    while stack:
        p = stack.pop()
        out.append(p)
        stack.extend(children.get(p, []))
    return out


def sample(pids):
    threads, cpu = 0, 0
    for p in pids:
        try:
            with open(f"/proc/{p}/stat") as f:
                fields = f.read().rsplit(")", 1)[1].split()
            cpu += int(fields[11]) + int(fields[12])  # utime + stime
            threads += int(fields[17])                 # num_threads
        except (OSError, IndexError, ValueError):
            pass
    return threads, cpu


def main():
    sep = sys.argv.index("--")
    log_path, cmd = sys.argv[1], sys.argv[sep + 1:]
    proc = subprocess.Popen(cmd)
    peak_threads, busy = 0, []
    last_cpu, last_t = 0, time.time()
    with open(log_path, "w") as log:
        log.write("elapsed_s\tthreads\tcores_busy\n")
        start = last_t
        while proc.poll() is None:
            time.sleep(INTERVAL)
            threads, cpu = sample(descendants(proc.pid))
            now = time.time()
            cores = max(0.0, (cpu - last_cpu) / TICKS / (now - last_t))
            last_cpu, last_t = cpu, now
            peak_threads = max(peak_threads, threads)
            busy.append(cores)
            log.write(f"{now - start:.0f}\t{threads}\t{cores:.2f}\n")
            log.flush()
    busy = busy[1:] or [0]
    print(f"MONITOR exit={proc.returncode} wall={time.time() - start:.0f}s "
          f"peak_threads={peak_threads} peak_cores={max(busy):.1f} "
          f"mean_cores={sum(busy) / len(busy):.1f} (machine has {os.cpu_count()})")
    sys.exit(proc.returncode)


if __name__ == "__main__":
    main()

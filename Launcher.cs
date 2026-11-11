using System;
using System.IO;
using System.Net;
using System.Drawing;
using System.Threading;
using System.Diagnostics;
using System.Windows.Forms;
using System.Security.Principal;

class Program {
    const string IMAGE_NAME = "hassancreator/biobio-app:1.0.0";
    const string IMAGE_DIGEST_FILE = "biobio-image.digest";
    const string CONTAINER_NAME = "biobio_desktop_instance";
    const string URL = "http://localhost:8501";

    static string DigestOnly(string value) {
        if (string.IsNullOrEmpty(value)) return "";
        int at = value.IndexOf('@');
        string digest = at >= 0 ? value.Substring(at + 1) : value;
        return digest.Trim();
    }

    // Fail closed when a release ships a digest sidecar: a modified or
    // re-tagged image must never run silently.
    static bool VerifyPinnedImage(string dockerExe, out string detail) {
        detail = "";
        string sidecar = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, IMAGE_DIGEST_FILE);
        if (!File.Exists(sidecar)) {
            detail = "development mode: no " + IMAGE_DIGEST_FILE;
            return true;
        }
        string expected = DigestOnly(File.ReadAllText(sidecar).Trim());
        if (expected.Length == 0) {
            detail = "digest sidecar is empty";
            return false;
        }
        string stdout, stderr;
        string args = "image inspect --format \"{{index .RepoDigests 0}}\" " + IMAGE_NAME;
        if (RunCmd(dockerExe, args, out stdout, out stderr) != 0) {
            detail = "docker image inspect failed: " + stderr;
            return false;
        }
        string actual = DigestOnly(stdout.Trim());
        if (actual != expected) {
            detail = "expected " + expected + " but found " + actual;
            return false;
        }
        detail = actual;
        return true;
    }

    static string GetDockerExe() {
        string defaultPath = @"C:\Program Files\Docker\Docker\resources\bin\docker.exe";
        if (File.Exists(defaultPath)) return defaultPath;
        return "docker";
    }

    static bool IsAdmin() {
        using (WindowsIdentity id = WindowsIdentity.GetCurrent()) {
            WindowsPrincipal pr = new WindowsPrincipal(id);
            return pr.IsInRole(WindowsBuiltInRole.Administrator);
        }
    }

    static int RunCmd(string file, string args, out string stdout, out string stderr) {
        stdout = "";
        stderr = "";
        try {
            ProcessStartInfo psi = new ProcessStartInfo(file, args);
            psi.RedirectStandardOutput = true;
            psi.RedirectStandardError = true;
            psi.UseShellExecute = false;
            psi.CreateNoWindow = true;
            using (Process p = Process.Start(psi)) {
                string outTmp = "";
                string errTmp = "";
                Thread tOut = new Thread(() => { outTmp = p.StandardOutput.ReadToEnd(); });
                Thread tErr = new Thread(() => { errTmp = p.StandardError.ReadToEnd(); });
                tOut.Start();
                tErr.Start();
                p.WaitForExit();
                tOut.Join();
                tErr.Join();
                stdout = outTmp;
                stderr = errTmp;
                return p.ExitCode;
            }
        } catch (Exception ex) {
            stderr = ex.Message;
            return -1;
        }
    }

    [STAThread]
    static void Main() {
        Application.EnableVisualStyles();
        Application.SetCompatibleTextRenderingDefault(false);

        string dockerExe = GetDockerExe();
        string outStr, errStr;
        int chk = RunCmd(dockerExe, "ps", out outStr, out errStr);

        if (chk != 0 && !IsAdmin()) {
            ProcessStartInfo elev = new ProcessStartInfo(Process.GetCurrentProcess().MainModule.FileName);
            elev.Verb = "runas";
            elev.UseShellExecute = true;
            try {
                Process.Start(elev);
                return;
            } catch {
                return;
            }
        }

        if (chk != 0) {
            MessageBox.Show(
                "يرجى فتح برنامج Docker Desktop أولاً والانتظار حتى يعمل المحرك (Engine running)، ثم أعد تشغيل البرنامج.\n\n" + errStr,
                "BIOBIO - تنبيه التشغيل",
                MessageBoxButtons.OK,
                MessageBoxIcon.Warning
            );
            return;
        }

        Form splash = new Form();
        splash.Text = "BIOBIO Suite";
        splash.Size = new Size(500, 180);
        splash.StartPosition = FormStartPosition.CenterScreen;
        splash.FormBorderStyle = FormBorderStyle.None;
        splash.BackColor = Color.FromArgb(15, 23, 42);

        Label titleLbl = new Label();
        titleLbl.Text = "BIOBIO — Molecular Docking & MD Suite";
        titleLbl.ForeColor = Color.FromArgb(56, 189, 248);
        titleLbl.Font = new Font("Segoe UI", 13, FontStyle.Bold);
        titleLbl.AutoSize = false;
        titleLbl.TextAlign = ContentAlignment.MiddleCenter;
        titleLbl.Dock = DockStyle.Top;
        titleLbl.Height = 70;

        Label statusLbl = new Label();
        statusLbl.Text = "جاري فحص وتجهيز المحرك الفيزيائي...";
        statusLbl.ForeColor = Color.White;
        statusLbl.Font = new Font("Segoe UI", 10, FontStyle.Regular);
        statusLbl.AutoSize = false;
        statusLbl.TextAlign = ContentAlignment.TopCenter;
        statusLbl.Dock = DockStyle.Fill;

        splash.Controls.Add(statusLbl);
        splash.Controls.Add(titleLbl);
        splash.Show();
        Application.DoEvents();

        if (RunCmd(dockerExe, "image inspect " + IMAGE_NAME, out outStr, out errStr) != 0) {
            statusLbl.Text = "جاري تحميل محرك BIOBIO لأول مرة من السحابة...\n(يرجى الانتظار ريثما يكتمل التحميل)";
            Application.DoEvents();
            int pullCode = -1;
            Thread pullThread = new Thread(() => {
                string po, pe;
                pullCode = RunCmd(dockerExe, "pull " + IMAGE_NAME, out po, out pe);
                errStr = pe;
            });
            pullThread.Start();
            while (pullThread.IsAlive) {
                Application.DoEvents();
                Thread.Sleep(200);
            }
            if (pullCode != 0) {
                splash.Close();
                MessageBox.Show("تعذر تحميل المحرك من السحابة:\n" + errStr, "BIOBIO Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }
        }

        statusLbl.Text = "جاري تشغيل المحرك الفيزيائي (GROMACS + Vina + ACPYPE)...";
        Application.DoEvents();

        string digestDetail;
        if (!VerifyPinnedImage(dockerExe, out digestDetail)) {
            splash.Close();
            MessageBox.Show(
                "فشل التحقق من بصمة صورة BIOBIO (digest).\n" +
                "هذا يعني أن الصورة تغيّرت أو أن الوسم أُعيد توجيهه، ولن يُسمح بالتشغيل.\n\n" +
                digestDetail,
                "BIOBIO - تحقق أمني",
                MessageBoxButtons.OK,
                MessageBoxIcon.Error
            );
            return;
        }

        string baseDir = AppDomain.CurrentDomain.BaseDirectory.TrimEnd('\\', '/');
        string dockingDir = Path.Combine(baseDir, "docking_runs");
        string mdDir = Path.Combine(baseDir, "md_runs");
        Directory.CreateDirectory(dockingDir);
        Directory.CreateDirectory(mdDir);

        RunCmd(dockerExe, "rm -f " + CONTAINER_NAME, out outStr, out errStr);

        string runArgs = string.Format(
            "run -d --name {0} -p 8501:8501 -v \"{1}:/app/docking_runs\" -v \"{2}:/app/md_runs\" {3}",
            CONTAINER_NAME, dockingDir, mdDir, IMAGE_NAME
        );

        if (RunCmd(dockerExe, runArgs, out outStr, out errStr) != 0) {
            splash.Close();
            MessageBox.Show("تعذر بدء حاوية البرنامج:\n" + errStr, "BIOBIO Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            return;
        }

        bool ready = false;
        for (int i = 0; i < 60; i++) {
            Application.DoEvents();
            try {
                HttpWebRequest req = (HttpWebRequest)WebRequest.Create(URL);
                req.Timeout = 1000;
                using (HttpWebResponse resp = (HttpWebResponse)req.GetResponse()) {
                    if (resp.StatusCode == HttpStatusCode.OK) {
                        ready = true;
                        break;
                    }
                }
            } catch {
                Thread.Sleep(500);
            }
        }

        splash.Hide();

        if (!ready) {
            RunCmd(dockerExe, "logs " + CONTAINER_NAME, out outStr, out errStr);
            MessageBox.Show("تأخر استجابة المحرك:\n" + outStr + "\n" + errStr, "BIOBIO Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            RunCmd(dockerExe, "rm -f " + CONTAINER_NAME, out outStr, out errStr);
            return;
        }

        // فتح نافذة البرنامج الأصلية عبر pywebview مباشرة بدون أي ملفات خارجية وبدون فتح متصفح Brave
        string pyCode = "import webview; webview.create_window('BIOBIO — Molecular Docking & MD Suite', 'http://localhost:8501', width=1440, height=920, min_size=(1024, 700)); webview.start()";
        ProcessStartInfo wpsi = new ProcessStartInfo("pythonw", "-c \"" + pyCode + "\"");
        wpsi.UseShellExecute = false;
        wpsi.CreateNoWindow = true;
        using (Process wp = Process.Start(wpsi)) {
            wp.WaitForExit();
        }

        RunCmd(dockerExe, "rm -f " + CONTAINER_NAME, out outStr, out errStr);
    }
}

import React, { useState } from "react";

function App() {
  const [smilesInput, setSmilesInput] = useState("CCO\nc1ccccc1\nCC(=O)Oc1ccccc1C(=O)O");
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const analyzeSmiles = async () => {
    setLoading(true);
    setError(null);
    setResults(null);

    const smilesArray = smilesInput.split("\n").map(s => s.trim()).filter(Boolean);

    try {
      const response = await fetch("http://127.0.0.1:8000/analyze", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles: smilesArray }),
      });

      if (!response.ok) {
        throw new Error(`API error: ${response.statusText}`);
      }

      const data = await response.json();
      setResults(data);
    } catch (err) {
      setError(err.message);
    }
    setLoading(false);
  };

  return (
    <div style={{ maxWidth: 600, margin: "auto", padding: 20, fontFamily: "Arial" }}>
      <h1>Chem-Cleanser Analyzer</h1>

      <textarea
        rows={6}
        value={smilesInput}
        onChange={(e) => setSmilesInput(e.target.value)}
        style={{ width: "100%", fontSize: 16, fontFamily: "monospace" }}
      />

      <button onClick={analyzeSmiles} disabled={loading} style={{ marginTop: 10, padding: "10px 20px" }}>
        {loading ? "Analyzing..." : "Analyze SMILES"}
      </button>

      {error && <p style={{ color: "red" }}>Error: {error}</p>}

      {results && (
        <div style={{ marginTop: 20 }}>
          <h2>Results</h2>
          <table border="1" cellPadding="6" style={{ borderCollapse: "collapse", width: "100%" }}>
            <thead>
              <tr>
                <th>SMILES</th>
                <th>Lipinski Pass</th>
                <th>Veber Pass</th>
                <th>MW Pass</th>
              </tr>
            </thead>
            <tbody>
              {results.smiles.map((smi, idx) => (
                <tr key={idx}>
                  <td>{smi}</td>
                  <td>{results.lipinski_pass[idx] ? "✅" : "❌"}</td>
                  <td>{results.veber_pass[idx] ? "✅" : "❌"}</td>
                  <td>{results.mw_pass ? (results.mw_pass[idx] ? "✅" : "❌") : "N/A"}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}

export default App;
